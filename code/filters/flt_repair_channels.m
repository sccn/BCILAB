function [signal,state] = flt_repair_channels(varargin)
% Repair (interpolate) broken channels online.
% [Signal,State] = flt_repair_channels(Signal,MinCorrelation,WindowLength,ProcessingDelay,ReferenceExtraction,CalibPrecision,LineNoiseAware,UseGPU)
%
% This is an automated artifact rejection function which temporarily interpolates channels based on
% the others when they become decorrelated (e.g. loose). The correlation is between the channel signal
% and the estimate of the channel signal according to other channels in a short moving window.
%
% In:
%   Signal : Continuous data set, assumed to be appropriately high-passed (e.g. >0.5Hz or with a 
%            0.5Hz - 2.0Hz transition band)
%
%   MinimumCorrelation  : A channel must have at least this high a correlation to its estimate based
%                         on other channels, otherwise it gets repaired. (default: 0.7)
%
%   WindowLength    : Window length over which correlation is calculated; this should be as short as
%                     possible while still being able to calculate accurate correlations (default: 0.5)
%
%   ProcessingDelay : Signal delay, in seconds. This can be anywhere between 0 and WindowLength/2;
%                     if this is too low some sharp-onset artifacts will leave brief glitches in the
%                     data. (default: 0.125).
%
%
%   The following are detail parameters that usually do not have to be tuned. If you can't get
%   the function to do what you want, you might consider adapting these to your data.
%
%   NumSamples : Number of RANSAC samples. This is the number of samples to generate in the random
%                sampling consensus process. (default: 20)
%
%   SubsetSize : Subset size. This is the size of the channel subsets to use, as number of channels
%                or a fraction of the total number of channels. (default: 0.25)
%
%   ReferenceExtraction : Group of sub-parameters to control how clean reference data shall be extracted 
%                         for calibration of the statistics (these are arguments to flt_clean_windows).
%                         Can also be set to 'off' to use the whole initialization data.
%
%       PowerTolerances: The minimum and maximum standard deviations within which the power of a channel
%                        must lie (relative to a robust estimate of the clean EEG power distribution in 
%                        the channel) for it to be considered "not bad". Default: [-3.5 5.5].
%
%       MaxBadChannels : The maximum number or fraction of bad channels that a retained window may still
%                        contain (more than this and it is removed). Reasonable range is 0.05 (very clean
%                        output) to 0.3 (very lax cleaning of only coarse artifacts). Default: 0.15.
%
%   InitializeOn : Initialize on time range. If a time range is given as [start,end], either in 
%                  seconds or as fractions of the whole data (if both <= 1), then where applicable 
%                  the filter will be initialized only on that subset of the whole data. As a 
%                  result, the filter will not have to be retrained in each cross-validation 
%                  fold. (default: [])
%
%   CalibPrecision : Block granularity for calibration. The data covariance will be estimated in
%                    blocks of this size and then robustly averaged (larger values admit smaller
%                    memory requirements). (default: 10)
%
%   LineNoiseAware : Whether the operation should be performed in a line-noise aware manner. If enabled,
%                    the correlation measure will not be affected by the presence or absence of line 
%                    noise. (default: true).
%
%   UseGPU : Whether to run on the GPU. Makes sense for offline processing if you have a GTX 780 or 
%            better. Default: false
%
%   State        :   previous filter state, as obtained by a previous execution of the filter on an
%                    immediately preceding data set (default: [])
%
% Out:
%   Signal : data set with bad channels removed
%
%   State  :  state of the filter, after it got applied
%
% Examples:
%   % use with defaults
%   eeg = flt_repair_channels(eeg);
%
%   % override the MinimumCorrelation default (making it more aggressive)
%   eeg = flt_clean_channels(eeg,0.7);
%
%   % override the MinimumCorrelation and the WindowLength, using name-value pairs
%   eeg = flt_clean_channels('Signal',eeg,'MinimumCorrelation',0.7, 'WindowLength',1.5);
%
% See also:
%   flt_clean_settings
%
% TODO:
%   The notch filter is a somewhat lacking solution.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-01-10

% flt_repair_channels_version<0.9.5> -- for the cache

if ~exp_beginfun('filter') return; end;

declare_properties('name','ChannelRepair', 'independent_channels',false, 'independent_trials','initialize_on');

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'min_corr','MinimumCorrelation'}, 0.7, [0 0.3 0.95 1], 'Minimum correlation between channels. This controls the aggressiveness of the filter: if the measure falls below this threshold for a channel in some time window, the window of that channel is considered bad.'), ...
    arg({'window_len','WindowLength'}, 0.5, [0 0.25 5 Inf], 'Window length to compute correlations. The length of the windows (in seconds) for which channel "badness" is computed, i.e. time granularity of the measure; ideally short enough to reasonably capture periods of artifacts, but no shorter (otherwise the statistic becomes too noisy).'),...
    arg({'processing_delay','ProcessingDelay'}, 0.05, [0 0 2.5 Inf], 'Signal delay, in seconds. This can be anywhere between 0 and WindowLength/2; if this is too low some sharp-onset artifacts will leave brief "knacks" in the data.'), ...
    arg({'num_samples','NumSamples'}, 20, uint32([1 10 100 10000]), 'Number of RANSAC samples. This is the number of samples to generate in the random sampling consensus process.'), ...
    arg({'subset_size','SubsetSize'}, 0.25, [0 Inf], 'Subset size. This is the size of the channel subsets to use, as number of channels or a fraction of the total number of channels.'), ...
    arg_subtoggle({'ref_extraction','ReferenceExtraction'},{'zthresholds',[-3.5 5.5]}, @flt_clean_windows, 'Reference data extraction. Group of sub-parameters to control how clean reference data shall be extracted for calibration of the statistics.'), ...
    arg({'initialize_on','InitializeOn'},[],[0 0 600 Inf],'Initialize on time range. If a time range is given as [start,end], either in seconds or as fractions of the whole data (if both <= 1), then where applicable the filter will be initialized only on that subset of the whole data. As a result, it will not have to be retrained in each cross-validation fold.','shape','row'),...
    arg({'calib_precision','CalibPrecision'}, 10, uint32([1 1 1000 10000]), 'Block granularity for calibration. The data covariance will be estimated in blocks of this size and then robustly averaged (larger values admit smaller memory requirements).','guru',true), ...
    arg({'linenoise_aware','LineNoiseAware'},true,[],'Line-noise aware processing. If enabled, a notch filter will be applied to ensure that line noise does not affect the correlation measure.','guru',true), ...
    arg({'chunk_length','ChunkLength'},50000,uint32([1 1000 100000 1000000000]), 'Maximum chunk length. Process the data in chunks of no larger than this (to avoid memory limitations).','guru',true), ...
    arg({'use_gpu','UseGPU'}, false, [], 'Whether to run on the GPU. Makes sense for offline processing if you have a GTX 780 or better.'), ...
    arg_deprecated({'history_fraction','HistoryFraction'}), arg_deprecated({'history_len','HistoryLength'}), arg_deprecated({'prefer_ica','PreferICAoverPCA'},false), ...
    arg_deprecated({'pca_flagquant','PCACleanliness'}),arg_deprecated({'pca_maxchannels','PCAForgiveChannels'}), ...
    arg_nogui({'state','State'}));

if subset_size < 1 %#ok<NODEF>
    subset_size = round(subset_size*size(signal.data,1)); end %#ok<NODEF>

signal.data = double(signal.data);
signal.data(~isfinite(signal.data(:))) = 0;

N = round(window_len*signal.srate);
P = round(processing_delay*signal.srate);

% initialize filter state if necessary
if isempty(state) %#ok<NODEF>
    if ~isempty(initialize_on)
        ref_section = exp_eval(set_selinterval(signal,initialize_on,quickif(all(initialize_on<=1),'fraction','seconds')));
    else
        ref_section = signal;
    end 
    if ref_extraction.arg_selection
        % extract a relatively clean section of data
        disp('Finding a clean section of the data...');
        ref_section = exp_eval_optimized(flt_clean_windows(ref_extraction,'signal',signal));
    end
    % determine the notch filter (if any)
    [state.B,state.A] = deal([]);
    if linenoise_aware
        try
            if signal.srate > 130
                [state.B,state.A] = yulewalk(8,[[0 45 50 55 60 65]*2/signal.srate 1],[1 1 0 1 0 1 1]);
            elseif signal.srate > 110
                [state.B,state.A] = yulewalk(4,[[0 45 50 55]*2/signal.srate 1],[1 1 0 1 1]);
            end
        catch %#ok<CTCH>
        end
    end
 
    % apply the notch filter and initialize the IIR filter state
    [C,S] = size(ref_section.data);
    if state.B
        [X,state.iir] = filter(state.B,state.A,ref_section.data,[],2); X = X';
    else
        X = ref_section.data';
    end
        
    % calculate the sample covariance matrices U (averaged in blocks of blocksize successive samples)
    if use_gpu
        try X = gpuArray(X); catch,end; end %#ok<CTCH>
    U = zeros(length(1:calib_precision:S),C*C);
    for k=1:calib_precision
        range = min(S,k:calib_precision:(S+k-1));
        U = U + reshape(bsxfun(@times,reshape(X(range,:),[],1,C),reshape(X(range,:),[],C,1)),size(U));
    end

    % get the reference mixing matrix M
    state.M = sqrtm(real(reshape(geometric_median(U/calib_precision),C,C)));
    
    % calculate randomly sampled reconstruction matrices
    for k=1:num_samples
        subset = randsample(1:C,subset_size);
        keep = false(C,1); keep(subset) = true;
        state.rand_samples{k} = real(state.M*pinv(bsxfun(@times,keep,state.M)))';
    end
    state.projector = horzcat(state.rand_samples{:});
    
    % initialize prior filter state by extrapolating available data into the past (if necessary)
    state.carry = repmat(2*signal.data(:,1),1,P) - signal.data(:,1+mod(((P+1):-1:2)-1,S));
    
    % and set other filters to zero
    [state.X_mov,state.Y_mov,state.XY_mov] = deal([]);
    % this is our breakage pattern database with associated reconstruction matrices
    state.patterns = [];
    state.reconstructors = {};
end

[C,S] = size(signal.data);

% pre-pend carry-over state
signal.data = [state.carry signal.data];

% split up the total sample range into k chunks that will fit in memory
splits = ceil(S/chunk_length);
if splits > 1
    fprintf('Now repairing channels in %i blocks',splits); end

for i=1:splits
    range = 1+floor((i-1)*S/splits) : min(S,floor(i*S/splits));
    if ~isempty(range)        
        % get new chunk of data (range shifted by processing_delay)
        X = signal.data(:,range+P);
        % get filtered data X
        if state.B
            [X,state.iir] = filter(state.B,state.A,X,state.iir,2); end
        % move it to the GPU if applicable
        if use_gpu && length(range) > 1000
            try X = gpuArray(X); catch,end; end %#ok<CTCH>
        
        % get the RANSAC'ed reconstruction Y
        Y = sort(reshape(X'*state.projector,size(X,2),C,[]),3);
        Y = Y(:,:,round(end/2))';
        
        % calc running variances E[X.^2] and E[Y.^2]
        [X_var,state.X_mov] = moving_average(N,X.^2,state.X_mov);
        [Y_var,state.Y_mov] = moving_average(N,Y.^2,state.Y_mov);
        % calc running covariance E[X.*Y]
        [XY_cov,state.XY_mov] = moving_average(N,X.*Y,state.XY_mov);
        % calc running correlation  (assuming zero means)
        XY_corr = XY_cov ./ (sqrt(abs(X_var)).*sqrt(abs(Y_var)));
        
        % determine matrix of bad channels
        bad_channels = XY_corr < min_corr;
        if use_gpu
            bad_channels = gather(bad_channels); end
            
        % get a mask of samples that need handling
        bad_samples = any(bad_channels);
        pattern_indices = zeros(size(bad_samples));
        % get the unique breakage patterns and a vector that encodes where the patterns fall
        [patterns,dummy,pattern_indices(bad_samples)] = unique(bad_channels(:,bad_samples)','rows'); %#ok<ASGLU>
        
        % for each such pattern...
        for p=1:size(patterns,1)
            patt = patterns(p,:);
            % generate the corresponding reconstruction matrix
            reconstruct = state.M*pinv(state.M(~patt,:));
            % apply it to the samples where the pattern occurs
            mask = p==pattern_indices;
            signal.data(patt,range(mask)) = reconstruct(patt,:) * signal.data(~patt,range(mask));
        end
    end
    if splits > 1
        fprintf('.'); end
end
if splits > 1
    fprintf('\n'); end

% annotate signal with filter delay
if ~isfield(signal.etc,'filter_delay')
    signal.etc.filter_delay = 0; end
signal.etc.filter_delay = signal.etc.filter_delay+round(processing_delay*signal.srate)/signal.srate;

% carry the look-ahead portion of the data over to the state (for successive calls)
state.carry = [state.carry signal.data(:,(end-P+1):end)];
state.carry = state.carry(:,(end-P+1):end);
if use_gpu
    state.iir = gather(state.iir);
    state.X_mov = gather(state.X_mov);
    state.Y_mov = gather(state.Y_mov);
    state.XY_mov = gather(state.XY_mov);
    state.projector = gather(state.projector);
    state.M = gather(state.M);
end

% finalize outputs
signal.data = signal.data(:,1:(end-P));

exp_endfun;



function [X,Zf] = moving_average(N,X,Zi)
% Run a moving-average filter along the second dimension of the data.
% [X,Zf] = moving_average(N,X,Zi)
%
% In:
%   N : filter length in samples
%   X : data matrix [#Channels x #Samples]
%   Zi : initial filter conditions (default: [])
%
% Out:
%   X : the filtered data
%   Zf : final filter conditions
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2012-01-10

if nargin <= 2 || isempty(Zi)
    Zi = zeros(size(X,1),N); end

% pre-pend initial state & get dimensions
Y = [Zi X]; M = size(Y,2);
% get alternating index vector (for additions & subtractions)
I = [1:M-N; 1+N:M];
% get sign vector (also alternating, and includes the scaling)
S = [-ones(1,M-N); ones(1,M-N)]/N;
% run moving average
X = cumsum(bsxfun(@times,Y(:,I(:)),S(:)'),2);
% read out result
X = X(:,2:2:end);

if nargout > 1
    Zf = [-(X(:,end)*N-Y(:,end-N+1)) Y(:,end-N+2:end)]; end



function Y = randsample(X,num)
Y = [];
while length(Y)<num
    pick = round(1 + (length(X)-1).*rand);
    Y(end+1) = X(pick);
    X(pick) = [];
end
