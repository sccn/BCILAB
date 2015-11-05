function [signal,state] = flt_repair_bursts(varargin)
% Repairs local peak artifacts using the ASR method.
% [Signal,State] = flt_repair_bursts(Signal,StandardDevCutoff,WindowLength,BlockSize,ProcessingDelay,MaxDimensions,SpectralParameters,CalibPrecision,UseGPU,ReferenceExtraction,State)
%
% This is an automated artifact rejection function that ensures that the data contains no events
% that have abnormally strong power; the subspaces on which those events occur are reconstructed 
% (interpolated) based on the rest of the EEG signal during these time periods.
%
% The basic principle is to first find a section of data that represents clean "reference" EEG and
% to compute statistics on there. Then, the function goes over the whole data in a sliding window
% and finds the subspaces in which there is activity that is more than a few standard deviations
% away from the reference EEG (this threshold is a tunable parameter). Once the function has found
% the bad subspaces it will treat them as missing data and reconstruct their content using a mixing
% matrix that was calculated on the clean data.
%
%
% In:
%   Signal : Continuous data set, assumed to be appropriately high-passed (e.g. >0.5Hz or with a 
%            0.5Hz - 2.0Hz transition band)
%
%   Cutoff : Standard deviation cutoff for removal of bursts (via ASR). Data portions whose variance
%            is larger than this threshold relative to the calibration data are considered missing
%            data and will be removed. The most aggressive value that can be used without losing too
%            much EEG is 2.5. A very conservative value would be 5. (default: 5)
%
%
%   The following are detail parameters that usually do not have to be tuned. If you cannot get
%   the function to do what you want, you might consider adapting these better to your data.
%
%   WindowLength : Length of the statistcs window, in seconds. This should not be much longer 
%                  than the time scale over which artifacts persist, but the number of samples in
%                  the window should not be smaller than 1.5x the number of channels. 
%                  (default: max(0.5,1.5*Signal.nbchan/Signal.srate))
%
%   StepSize : Step size for processing. The reprojection matrix will be updated every this many
%              samples and a blended matrix is used for the in-between samples. If empty this will
%              be set the WindowLength/2 in samples. (default: 1/3)
%
%   ProcessingDelay : Signal delay, in seconds. This can be anywhere between 0 and WindowLength/2;
%                     if this is too low some sharp-onset artifacts will leave brief "knacks" in the
%                     data. (default: 0.125).
%
%   MaxDimensions : Maximum dimensionality to reconstruct. Up to this many dimensions (or up to this 
%                   fraction of dimensions) can be reconstructed for a given data segment. This is
%                   since the lower eigenvalues are usually not estimated very well. Default: 2/3.
%
%   SpectralParameters: Parameters for the spectrum-shaping filter that reweights frequencies
%                       according to their relevance for artifact detection. This is a cell
%                       array of {Frequencies,Amplitudes,Order} for a Yule-Walker IIR filter
%                       design. Note: These settings may give you a different frequency response from
%                       what you expect, so check what you generate. Also, the default value is not 
%                       recommended for >512Hz sampling rate. 
%                       (default: {[0 2 3 13 16 40 80],[3 0.75 0.33 0.33 1 1 3],8})
%
%   CalibPrecision : Block size for calculating the robust data covariance, in samples; allows to 
%                    reduce the memory requirements of the robust estimator by this factor (down to
%                    Channels x Channels x Samples x 16 / Blocksize bytes). The value is further
%                    incremented if necessary to fit the data in memory. (default: 10).
%
%   UseGPU : Whether to run on the GPU. Makes sense for offline processing if you have a GTX 780 or 
%            better. Default: false
%
%   InitializeOn : Initialize on time range. If a time range is given as [start,end], either in 
%                  seconds or as fractions of the whole data (if both <= 1), then where applicable 
%                  the filter will be initialized only on that subset of the whole data. As a 
%                  result, the filter will not have to be retrained in each cross-validation 
%                  fold. (default: [])
%
%   ReferenceExtraction : Group of sub-parameters to control how clean reference data shall be extracted 
%                         for calibration of the statistics (these are arguments to flt_clean_windows).
%                         Can also be set to 'off' to select the entire initialization data.
%
%       PowerTolerances: The minimum and maximum standard deviations within which the power of a channel
%                        must lie (relative to a robust estimate of the clean EEG power distribution in 
%                        the channel) for it to be considered "not bad". Default: [-3.5 5.5].
%
%       MaxBadChannels : The maximum number or fraction of bad channels that a retained window may still
%                        contain (more than this and it is removed). Reasonable range is 0.05 (very clean
%                        output) to 0.3 (very lax cleaning of only coarse artifacts). Default: 0.15.
%
%   State : previous filter state, as obtained by a previous execution of the filter on an immediately 
%           preceding data set. If this is a string of the form ".fieldname" for a valid field name
%           in the input data set, the function will take this as its initial state and not
%           calibrate again. This is useful for carrying resting-state calibration over into a
%           regular BCILAB analysis, e.g. a cross-validation. (default: [])
%
% Out:
%   Signal : data set with local peaks removed
%
%   State  :  state of the filter, after it got applied
%
% Notes:
%   When used online, this function gives slightly different results from offline processing since
%   parameters are also updated at the beginning and end of each chunk; therefore, transitions can be 
%   as short as the chunk length. The behavior should nevertheless be qualitatively equivalent to
%   offline processing.
%
% See also:
%   flt_clean_settings
% 
% Examples:
%   % use the defaults
%   eeg = flt_repair_bursts(eeg);
%
%   % use a 1-second window length and 0.5 seconds delay / look-ahead (e.g. for very high channel counts)
%   eeg = flt_repair_bursts('Signal',eeg, 'WindowLength',1, 'ProcessingDelay',0.5)
%
%   % use some very aggressive settings (not generally advisable)
%   eeg = flt_repair_bursts('Signal',eeg, 'StandardDevCutoff',2.5,'MaxDimensions',0.75);
%
%   % use a very conservative setting
%   eeg = flt_repair_bursts(eeg,4);
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-08-25

% flt_repair_bursts _version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

declare_properties('name','BurstRepair', 'independent_channels',false, 'independent_trials','initialize_on');

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'stddev_cutoff','Cutoff','StandardDevCutoff','cutoff'}, 5, [0 2 15 Inf], 'StdDev cutoff for rejection. Data segments whose variance is beyond this cutoff from the distribution of variance across the recording are considered missing data.'), ...
    arg({'window_len','WindowLength'}, 0.5, [0 0.125 2 Inf], 'Window length to compute signal power. In seconds, ideally short enough to reasonably isolate artifacts, but no shorter, to keep statistics healthy.'), ...
    arg({'block_size','BlockSize'}, 1/3, [], 'Block granularity for processing. The reprojection matrix will be updated every this many samples. If this is below 1, it is assumed to be in seconds.','guru',true), ...
    arg({'processing_delay','ProcessingDelay'}, 0.125, [0 Inf], 'Signal delay, in seconds. This can be anywhere between 0 and WindowLength/2; if this is too low some sharp-onset artifacts will leave brief "knacks" in the data.'), ...
    arg({'max_dimensions','MaxDimensions'}, 0.66, [0 Inf], 'Maximum dimensionality to reconstruct. Up to this many dimensions (or up to this fraction of dimensions) can be reconstructed for a given data segment. This is useful because the lower eigenvalues are usually not estimated very well.','guru',true), ...
    arg({'spectral_params','SpectralParameters'},{[0 2 3 13 16 40 80],[3 0.75 0.33 0.33 1 1 3],8},[],'Spectral weighting parameters. Parameters for the spectrum-shaping filter that reweights frequencies according to their relevance for artifact detection. This is a cell array of {Frequencies,Amplitudes,Order} for a Yule-Walker IIR filter design. Note: These settings are very sensitive, do not change unless you know exactly what you''re doing; specifically, the default it not recommended for >512 Hz srate.','guru',true), ...
    arg({'calib_precision','CalibPrecision'}, 10, uint32([1 1 1000 10000]), 'Block granularity for calibration. The data covariance will be estimated in blocks of this size and then robustly averaged (larger values admit smaller memory requirements).','guru',true), ...
    arg({'use_gpu','UseGPU'}, false, [], 'Whether to run on the GPU. Makes sense for offline processing if you have a card with double-precision performance of a GTX Titan or better.'), ...
    arg({'processing_mode','ProcessingMode'}, 'remove', {'remove','keep','separate'}, 'Processing mode. Determines what to do after artifacts have been identified.'), ...
    arg({'separate_noise','SeparateNoise'}, 0.001, [0, 0.00001, 0.01, 1], 'Noise level for separated artifacts. This is the relative amplitude of random noise mixed into the (low-rank) artifact signal, if ProcessingMode is separate, to keep downstream algorithms from yielding nonfinite numbers.'), ...
    arg({'initialize_on','InitializeOn'},[],[0 0 600 Inf],'Initialize on time range. If a time range is given as [start,end], either in seconds or as fractions of the whole data (if both <= 1), then where applicable the filter will be initialized only on that subset of the whole data. As a result, it will not have to be retrained in each cross-validation fold.','shape','row'),...
    arg_subtoggle({'ref_extraction','ReferenceExtraction'},{'zthresholds',[-3.5 5.5]}, @flt_clean_windows, 'Reference data extraction. Group of sub-parameters to control how clean reference data shall be extracted for calibration of the statistics.'), ...
    arg_deprecated({'ref_wndlen','ReferenceWindowLength'},[],[],'Selection granularity for reference data. When clean data is extracted from the recording for calibrating the cutoff thresholds the data is chopped up into windows of this length. This parameters has been deprecated in favor of ReferenceExtraction.WindowLength.'),...
    arg_deprecated({'ref_maxbadchannels','ReferenceMaxBadChannels'}, [], [], 'Maximum bad-channel fraction on reference data. The maximum fraction of bad channels per time window of the data that is used as clean reference EEG on which statistics are based. Instead of a number one may also directly pass in a data set that contains clean reference data (for example a minute of resting EEG). A lower values makes this criterion more aggressive. Reasonable range: 0.05 (very aggressive) to 0.3 (quite lax). This parameters has been deprecated in favor of ReferenceExtraction.MaxBadChannels.'), ...
    arg_deprecated({'ref_tolerances','ReferenceTolerances'},[],[],'Power tolerances on reference data. These are the power tolerances beyond which a channel in the clean reference data is considered "bad", in standard deviations relative to a robust EEG power distribution (lower and upper bound). This parameters has been deprecated in favor of ReferenceExtraction.PowerTolerances.'), ...
    arg_nogui({'state','State'}));

% handle legacy parameters
if ~isempty(ref_wndlen)
    ref_extraction.window_length = ref_wndlen; end
if ~isempty(ref_maxbadchannels)
    ref_extraction.max_bad_channels = ref_maxbadchannels; end
if ~isempty(ref_tolerances)
    ref_extraction.zthresholds = ref_tolerances; end

% ensure that the window length is large enough
window_len = max(window_len,1.5*signal.nbchan/signal.srate); %#ok<*NODEF>
% ensure that the spectral filter does not exceed the Nyquist rate of the signal
spectral_params{1} = min(spectral_params{1},signal.srate/2-1);
% set a reasonable processing delay/look-ahead
if isempty(processing_delay)
    processing_delay = 0.125; end
% also set a reasonable block size...
if isempty(block_size)
    block_size = 1/3; end
block_size = round(signal.srate*block_size);

% check if there is a calibration model already in a given field of the data
if ~isempty(state) && ischar(state)
    if isfield(signal,state(2:end))
        state = signal.(state(2:end));
        disp('Using cleaning model already contained in the data.');
    else
        state = []; % otherwise disregard the field
    end
end
    
% initialize filter state if necessary
if isempty(state)
    % design a spectral filter
    try
        [state.shaping{1},state.shaping{2}] = hlp_diskcache('filterdesign',@yulewalk,spectral_params{3},[2*spectral_params{1}/signal.srate 1],spectral_params{2}([1:end end]));
    catch e
        if strcmp(e.identifier,'MATLAB:UndefinedFunction')
            disp('The yulewalk function was not found to design a spectral filter; trying to fall back to pre-computed filter.');            
            state.shaping = {[],[]};
        else
            rethrow(e);
        end
    end
    if ~isempty(initialize_on)
        ref_section = exp_eval(set_selinterval(signal,initialize_on,quickif(all(initialize_on<=1),'fraction','seconds')));
    else
        ref_section = signal;
    end 
    if ref_extraction.arg_selection
        % extract reference data subset without gross artifacts
        disp('Finding a clean section of the data...');
        ref_section = exp_eval(flt_clean_windows('signal',ref_section,ref_extraction));
    end
    % calibrate on it
    disp('Estimating statistics...');
    state.asr = hlp_diskcache('filterdesign',@asr_calibrate,ref_section.data,ref_section.srate,stddev_cutoff,calib_precision,state.shaping{:}); clear ref_section;
end


if ~strcmp(processing_mode,'remove')
    old_data = signal.data;
    old_carry = state.asr.carry;
    if isempty(old_carry)
        P = round(processing_delay*signal.srate);
        S = size(old_data,2);
        old_carry = repmat(2*old_data(:,1),1,P) - old_data(:,1+mod(((P+1):-1:2)-1,S)); 
    end
end

% do the actual processing
[signal.data,state.asr] = asr_process(signal.data,signal.srate,state.asr,window_len,processing_delay,block_size,max_dimensions,[],use_gpu);

if strcmp(processing_mode,'remove')
    % nothing to do
else
    old_data = [old_carry old_data(:,1:end-size(state.asr.carry,2))];
    if strcmp(processing_mode,'keep')
        % keep the artifacts (difference between original and cleaned signal), remove the EEG
        signal.data = signal.data - old_data;
    elseif strcmp(processing_mode,'separate')
        % separate artifacts from EEG, keep both
        data_std = median(std(signal.data,0,2));
        signal.data = [signal.data; (signal.data - old_data) + data_std*separate_noise*randn(size(signal.data))];
        signal.chanlocs = [signal.chanlocs signal.chanlocs];
        for c=signal.nbchan:length(signal.chanlocs)
            signal.chanlocs(c).labels = [signal.chanlocs(c).labels '_art']; 
            signal.chanlocs(c).type = [signal.chanlocs(c).type '_art']; 
        end        
        signal.nbchan = length(signal.chanlocs);
    else
        error('Unsupported processing mode: %s',hlp_tostring(processing_mode));
    end
end

if ~isfield(signal.etc,'filter_delay')
    signal.etc.filter_delay = 0; end
signal.etc.filter_delay = signal.etc.filter_delay+round(processing_delay*signal.srate)/signal.srate;

exp_endfun;
