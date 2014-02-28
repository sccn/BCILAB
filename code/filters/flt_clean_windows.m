function [signal,sample_mask] = flt_clean_windows(varargin)
% Remove periods with abnormally high/low-amplitude content from continuous data.
% [Signal,Mask] = flt_clean_windows(Signal,PowerTolerances,MaxBadChannels,WindowLength,ParameterFitting,KeepMetadata,FlaggedQuantile,MinAffectedChannels,MaxIgnoredChannels)
%
% This function cuts segments from the data which contain high/low-amplitude artifacts.
% Specifically, any windows with more than a certain fraction of "bad" channels are removed, where a
% channel is bad in a given window if its amplitude in the window is above or below a given
% upper/lower threshold (in standard deviations from a robust estimate of the EEG amplitude
% distribution for the channel).
%
% The function also supports several legacy methods (for backwards compatibility).
%
% In:
%   Signal         : Continuous data set, assumed to be appropriately high-passed (e.g. >1Hz or
%                    0.5Hz - 1.0Hz transition band)
%
%   PowerTolerances: The minimum and maximum standard deviations within which the power of a channel
%                    must lie (relative to a robust estimate of the clean EEG power distribution in 
%                    the channel) for it to be considered "not bad". Default: [-3.5 5].
%
%   MaxBadChannels : The maximum number or fraction of bad channels that a retained window may still
%                    contain (more than this and it is removed). Reasonable range is 0.05 (very clean
%                    output) to 0.3 (very lax cleaning of only coarse artifacts). Default: 0.15.
%
%
%
%   The following are detail parameters that usually do not have to be tuned. If you can't get
%   the function to do what you want, you might consider adapting these to your data.
%
%   WindowLength : Window length that is used to check the data for artifact content. This is 
%                  ideally as long as the expected time scale of the artifacts but short enough to 
%				   allow for several 1000 windows to compute statistics over. Default: 0.5.
%
%   WindowOverlap : Window overlap fraction. The fraction of two successive windows that overlaps.
%                   Higher overlap ensures that fewer artifact portions are going to be missed (but
%                   is slower). (default: 0.66)
%
%   ParameterFitting : Group of sub-arguments that govern how EEG distribution parameters should be fit 
%
%       MaxDropoutFraction : Maximum fraction of values in X that can be subject to
%                            signal dropouts (e.g., sensor unplugged) (default: 0.1)
%
%       MinCleanFraction : Minimum fraction of values in X that needs to be clean
%                          (default: 0.25)
%
%       FitQuantiles : Quantile range [upper,lower] of the truncated Gaussian distribution 
%                      that shall be fit to the EEG contents (default: [0.022 0.6])
%
%       StepSizes : Step size of the grid search, in quantiles; separately for [lower,upper] edge of the
%                   truncated Gaussian (default: [0.001 0.01])
%
%       ShapeRange : Shape parameter range. Search range for the shape parameter of the generalized
%                    Gaussian distribution used to fit clean EEG. (default: 1.5:0.1:3.5)
%
%`  KeepMetadata    : boolean; whether meta data of EEG struct (such as events, ICA decomposition
%                     etc.) should be returned. If true, meta data is returned. Returning meta data
%                     is quite slow. (default: true)
%
%
%
%	The following are legacy parameters to enable previous behaviors.
%
%   FlaggedQuantile : upper quantile of the per-channel windows that should be flagged for potential
%                     removal (removed if flagged in all except for some possibly bad channels);
%                     controls the aggressiveness of the rejection; if two numbers are specified,
%                     the first is the lower quantile and the second is the upper quantile to be
%                     removed (default: [], was: 0.15)
%
%   MinAffectedChannels : if for a time window more than this number (or ratio) of channels are
%                         affected (i.e. flagged), the window will be considered "bad". (default: [], was: 0.5)
%
%   MaxIgnoredChannels : maximum number (or ratio) of channels which can deliver arbitrary data
%                        without affecting the outcome. (default: [], was: 0.3)
%
% Out:
%   Signal : data set with bad time periods removed.
%
%   Mask   : mask of retained samples (logical array)
%
% Examples:
%   % use the defaults
%   eeg = flt_clean_windows(eeg);
%
%   % use a more aggressive threshold and a custom window length
%   eeg = flt_clean_windows(eeg,0.25,0.5);
%
%   % use the default, but keep the meta-data (i.e. events, etc) - also, pass all parameters by name
%   eeg = flt_clean_windows('Signal',eeg,'KeepMetadata',true);
%
% See also:
%   flt_clean_channels, flt_clean_settings
% 
% Notes:
%   This method has no effect when used online.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-07-06

% flt_clean_windows_version<1.02> -- for the cache

if ~exp_beginfun('editing') return; end;

declare_properties('name','WindowCleaning', 'independent_channels',false, 'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'zthresholds','PowerTolerances'}, [-3.5 5], [], 'Z score thresholds. The minimum and maximum standard deviations within which the power of a channel must lie (relative to a robust estimate of the clean EEG power distribution in the channel) for it to be considered "not bad".'), ...
    arg({'max_bad_channels','MaxBadChannels'}, 0.15, [], 'Maximum fraction of bad channels. The maximum number or fraction of bad channels that a retained window may still contain (more than this and it is removed). Reasonable range is 0.05 (very clean output) to 0.3 (very lax cleaning of only coarse artifacts).'), ...
    arg({'window_len','WindowLength','window_length'}, 0.5, [], 'Window length. Window length that is used to check the data for artifact content. This is ideally as long as the expected time scale of the artifacts but short enough to allow for several 1000 windows to compute statistics over.'), ...
    arg({'window_overlap','WindowOverlap'}, 0.66, [0 1], 'Window overlap fraction. The fraction of two successive windows that overlaps. Higher overlap ensures that fewer artifact portions are going to be missed (but is slower).'), ...
    arg_sub({'fit_params','ParameterFitting','parameter_fitting'}, {}, { ...
        arg({'max_dropout_fraction','MaxDropoutFraction'}, 0.1, [0 1], 'Maximum fraction that can have dropouts. This is the maximum fraction of time windows that may have arbitrarily low amplitude (e.g., due to the sensors being unplugged).'), ...
        arg({'min_clean_fraction','MinCleanFraction'}, 0.25, [0 1], 'Minimum fraction that needs to be clean. This is the minimum fraction of time windows that need to contain essentially uncontaminated EEG.'), ...
        arg({'truncate_quant','TruncateQuantile'}, [0.022 0.6], [0 1], 'Truncated Gaussian quantile. Quantile range [upper,lower] of the truncated Gaussian distribution that shall be fit to the EEG contents.','guru',true), ...
        arg({'step_sizes','StepSizes'}, [0.01 0.01], [0.00001 0.2], 'Grid search stepping. Step size of the grid search, in quantiles; separately for [lower,upper] edge of the truncated Gaussian. The lower edge has finer stepping because the clean data density is assumed to be lower there, so small changes in quantile amount to large changes in data space.','guru',true), ...
        arg({'shape_range','ShapeRange'}, 1.7:0.15:3.5, [], 'Shape parameter range. Search range for the shape parameter of the generalized Gaussian distribution used to fit clean EEG.','guru',true), ...
        arg_deprecated({'num_bins','NumBins'},50,[],'This parameter is now auto-determined.') ...
        arg_deprecated({'fit_quantiles','FitQuantiles'}, [0.022 0.6], [], 'Old Truncated Gaussian quantile. This was allowed to be a matrix holding an ensemble of alternative quantiles -- that is no longer needed/supported (instead the replacement parameter TruncateQuantile should be used).','shape','matrix'), ...
    }, 'Parameter fitting details. Group of sub-arguments that govern how EEG distribution parameters should be fit.'), ...
    arg({'keep_metadata','KeepMetadata'}, true, [], 'Retain metadata of EEG set. Retaining meta data (events, ICA decomposition, etc.) is quite slow.'), ...
    arg({'debugdisplay','DebugDisplay'}, false, [], 'Enable debug display. Plots the histogram and fit for each channel.'), ...
    arg_deprecated({'flag_quantile','FlaggedQuantile'}, [], [], 'Legacy parameter for pre-2012 methods. Quantile of data windows flagged for removal. Windows are emoved if flagged in all except for some possibly bad channels, controls the aggressiveness of the rejection.'), ...
    arg_deprecated({'min_badchans','MinAffectedChannels'}, [], [], 'Legacy parameter for 2011 method. This is the minimum number of channels that need to be affected for the window to be considered "bad".'), ...
    arg_deprecated({'ignored_chans','MaxIgnoredChannels'}, [], [], 'Legacy parameter for 2010 method. Maximum number or ratio of channels to ignore. These can contain arbitrary data without affecting the outcome (e.g., can be unplugged channels).'));

if ~isempty(max_bad_channels) && max_bad_channels > 0 && max_bad_channels < 1 %#ok<*NODEF>
    max_bad_channels = round(size(signal.data,1)*max_bad_channels); end
if ~isempty(flag_quantile) && isscalar(flag_quantile)
    flag_quantile = [0 flag_quantile]; end
if ~isempty(min_badchans) && min_badchans > 0 && min_badchans < 1 %#ok<*NODEF>
    min_badchans = round(size(signal.data,1)*min_badchans); end
if ~isempty(ignored_chans) && ignored_chans > 0 && ignored_chans < 1
    ignored_chans = round(size(signal.data,1)*ignored_chans); end

signal.data = double(signal.data);
[C,S] = size(signal.data);
N = window_len*signal.srate;
wnd = 0:N-1;
offsets = round(1:N*(1-window_overlap):S-N);
W = length(offsets);

if max_bad_channels < 1
    max_bad_channels = round(max_bad_channels*C); end
if max_bad_channels >= C
    sample_mask = true(1,S); return; end

fprintf('Determining time window rejection thresholds...');
if debugdisplay
    figure; end
% for each channel...
for c = C:-1:1
    % compute RMS amplitude for each window...
    X = signal.data(c,:).^2;
    X = sqrt(sum(X(bsxfun(@plus,offsets,wnd')))/N);
    % robustly fit a distribution to the clean EEG part
    if isempty(flag_quantile)
        [mu,sig,alpha,beta] = hlp_diskcache('filterdesign',@fit_eeg_distribution,X, ...
            fit_params.min_clean_fraction, fit_params.max_dropout_fraction, ...
            fit_params.truncate_quant, fit_params.step_sizes,fit_params.shape_range); %#ok<NASGU,ASGLU>
        if debugdisplay
            clf; 
            show_fit_quality(X,mu,sig,alpha,beta); 
            title(sprintf('%i@(%.2f,%.2f,%.2f,%.2f)',c,mu,sig,alpha,beta)); 
            drawnow; 
            waitforbuttonpress;
        end
    else
        [mu,sig] = deal(0,1);
    end
    % calculate z scores relative to that
    wz(c,:) = (X - mu)/sig;
end
disp('done.');

% determine which windows to retain
if isempty(flag_quantile)
    % best known method (2013): based on z thresholds

    % sort z scores into quantiles
    swz = sort(wz);
    % determine which windows to remove
    remove_mask = false(1,size(swz,2));
    if max(zthresholds)>0
        remove_mask(swz(end-max_bad_channels,:) > max(zthresholds)) = true; end
    if min(zthresholds)<0
        remove_mask(swz(1+max_bad_channels,:) < min(zthresholds)) = true; end
    removed_windows = find(remove_mask);

else
    % legacy methods for backwards compatibility
    [dummy,i] = sort(wz'); %#ok<TRSRT,ASGLU>
    % find retained window indices per channel
    retained_quantiles = i(1+floor(W*flag_quantile(1)):round(W*(1-flag_quantile(2))),:)';
    % flag them in a Channels x Windows matrix (this can be neatly visualized)
    retain_mask = zeros(C,W);
    for c = 1:C
        retain_mask(c,retained_quantiles(c,:)) = 1; end
    % find retained windows
    if ~isempty(min_badchans) 
        % previous method (2011): based on quantiles
        removed_windows = find(sum(1-retain_mask) > min_badchans);
    elseif ~isempty(ignored_chans)
        % oldest method (2010): attempt at being robust against broken channels
        removed_windows = find(sum(retain_mask) <= ignored_chans);
    else
        error('BCILAB:flt_clean_windows:legacy_options','By setting the flag_quantile parameter you switch to legacy methods; these require that you also set the min_badchans (former default: 0.5) or ignored_chans (former default: 0.3) parameter.');
    end
end


% find indices of samples to remove
removed_samples = repmat(offsets(removed_windows)',1,length(wnd))+repmat(wnd,length(removed_windows),1);
% mask them out
sample_mask = true(1,S); 
sample_mask(removed_samples(:)) = false;
fprintf('Keeping %.1f%% (%.0f seconds) of the data.\n',100*(mean(sample_mask)),nnz(sample_mask)/signal.srate);

if keep_metadata
    % retain the masked data, update meta-data appropriately
    signal = exp_eval(set_selinterval(signal,sample_mask,'range'));
    
    if isfield(signal.etc,'epoch_bounds') && isfield(signal.event,'target')
        targets = find(~cellfun('isempty',{signal.event.target}));
        retain = targets;
        % further restrict the set of retained events: generate epoch index range, in samples
        eporange = round(signal.etc.epoch_bounds(1)*signal.srate) : round(signal.etc.epoch_bounds(2)*signal.srate);
        
        if ~isempty(eporange)
            % prune events that exceed the data set boundaries
            lats = round([signal.event(retain).latency]);
            retain(lats+eporange(1)<1 | lats+eporange(end)>signal.pnts) = [];
            
            % generate a sparse mask of boundary events
            boundlats = min(signal.pnts,max(1,round([signal.event(strcmp({signal.event.type},'boundary')).latency])));
            if ~isempty(boundlats)
                boundmask = sparse(ones(1,length(boundlats)),boundlats,1,1,signal.pnts);
                
                % prune events that intersect the boundary mask
                lats = round([signal.event(retain).latency]);
                if ~isempty(lats)
                    retain(any(boundmask(bsxfun(@plus,eporange',lats)))) = []; end
                
                % now remove them
                remove = setdiff(targets,retain);
                signal.event(remove) = [];
            end
        end
    end
    
else
    % retain the masked data, clear all events or other aggregated data
    signal = exp_eval(set_new(signal,'data',signal.data(:,sample_mask),'icaact',[],'event',signal.event([]),'urevent',signal.urevent([]), ...
        'epoch',[],'reject',[],'stats',[],'specdata',[],'specicaact',[]));
end

% annotate the data with the mask of what was removed (e.g., for visualization)
if isfield(signal.etc,'clean_sample_mask') && nnz(signal.etc.clean_sample_mask) == S
    signal.etc.clean_sample_mask(signal.etc.clean_sample_mask) = sample_mask;
else
    signal.etc.clean_sample_mask = sample_mask;
end

exp_endfun;