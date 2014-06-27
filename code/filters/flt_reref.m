function signal = flt_reref(varargin)
% Re-references the data to a new (set of) channel(s) or the average of all channels.
% Signal = flt_reref(Signal, ReferenceChannels, ExcludeChannels, KeepReference, ReferenceType, HuberCutoff, HuberIterations)
%
% Re-referencing is a spatial filter in which the (instantaneous) mean signal of some "reference"
% channels is subtracted from the signal of each channel. This allows to dampen externally-induced
% signals (e.g. electromagnetic interference) which are captured by all sensors, including those
% that are references. A frequently used variant is the average reference, in which the mean of all
% channels is subtracted from every channel. Dedicated reference channels are usually behind the
% ears ("mastoids"), on the nose bone ("nasion") or sometimes on the scalp, whereas this cancels
% brain signals out of nearby channels. Many BCI paradigms do not require an explicit
% re-referencing, since their spatial filters are adaptively optimized and effectively incorporate
% re-referencing to the degree that it is necessary.
%
% This function also supports robust re-referencing, where artifacts in few channels will not affect
% the outcome disproportionally.
%
% In:
%   Signal            :  Epoched or continuous data set.
%
%   ReferenceChannels :  Cell array of reference channels to which the signal data should be
%                        referenced, or [] for average reference. (default: [])
%
%   ExcludeChannels   : Cell array of channels to exclude from the re-referencing computation
%                       (default: [])
%
%   KeepReference     : Whether to keep the reference channel (default: false)
%                       If true and multiple channels were selected, a new channel named CAR is
%                       appended.
%
%   ReferenceType     : Type of referencing to apply, can be one of the following:
%                       * 'mean': subtract average of reference channels (default)
%                       * 'median': subtract median of reference channels
%                       * 'huber': subtract robust mean (under the Huber loss) of reference channels;
%
%   HuberCutoff       : cutoff/transition parameter for huber loss; if [], this is set
%                       to one (robust) standard deviation of the signal
%
%   HuberIterations   : number of iterations to compute the Huber fit; a larger number can tolerate
%                       larger outliers (default: 100)
%
% Out:
%   Signal : Re-referenced data set.
%
% Examples:
%   % do a common average reference
%   eeg = flt_reref(eeg)
%
%   % re-reference to the TP7 and TP8 electrodes
%   eeg = flt_reref(eeg,{'TP7','TP8'})
%
%   % do a common average reference, but exclude some assumed 'EOGV' and 'EOGH' channels
%   eeg = flt_reref(eeg,[],{'EOGV','EOGH'})
%
%   % as previous, but passing the arguments by name
%   eeg = flt_reref('Signal',eeg,'ExcludeChannels',{'EOGV','EOGH'})
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% flt_reref_version<1.01> -- for the cache

if ~exp_beginfun('filter') return; end

declare_properties('name',{'Rereferencing','ref','Rereference'}, 'independent_channels',false,'independent_trials',true);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'ref_chn','ReferenceChannels','chn'}, [], [], 'Cell array of reference channels. The signal data is be referenced to these, defaults to average reference if empty.','type','cellstr','shape','row'),...
    arg({'exclude_chn','ExcludeChannels','exclude'}, [], [], 'Cell array of channels to exclude.','type','cellstr','shape','row'),...
    arg({'keepref','KeepReference'}, false, [], 'Keep the reference channel.'), ...
    arg({'ref_type','ReferenceType'}, 'mean', {'mean','median','huber'}, 'Type of reference. The traditional average reference operation uses the mean. If this is set to median, the median of the reference channels will be removed, and if set to huber the robust mean under the Huber loss will be removed (slower than median but statistically a more efficient estimator).'), ...
    arg({'huber_cut','HuberCutoff'}, [], [], 'Cutoff for huber function. If left empty this is set to one (robust) standard deviation of the signal.','shape','scalar','guru',true), ...
    arg({'huber_iters','HuberIterations'}, 100, uint32([1 10 1000 10000]), 'Iterations for huber fitting. A larger number yields tolerance to larger outliers.','guru',true));

if ~isempty(ref_chn)
    ref_chns = set_chanid(signal,ref_chn);
else
    ref_chns = 1:size(signal.data,1);
end
    
if ~isempty(exclude_chn)
    ref_chns = setdiff(ref_chns,set_chanid(signal,exclude_chn)); end

switch ref_type
    case 'mean'
        ref_signal = mean(signal.data(ref_chns,:));
    case 'median'
        ref_signal = median(signal.data(ref_chns,:));
    case 'huber'
        if isempty(huber_cut)
            huber_cut = median(median(abs(bsxfun(@minus,signal.data,median(signal.data,2))),2))*1.4826; end
        ref_signal = robust_mean(signal.data(ref_chns,:)/huber_cut,1,huber_iters)*huber_cut;
    otherwise
        error('Unsupported reference type.');
end
    
signal.data = bsxfun(@minus,signal.data,ref_signal);

if ~keepref && ~isempty(ref_chn)
    retain_channels = true(1,size(signal.data,1)); 
    retain_channels(ref_chns) = false;
    signal.data = signal.data(retain_channels,:,:,:,:,:,:,:);
    signal.chanlocs = signal.chanlocs(retain_channels);
    signal.nbchan = size(signal.data,1);    
end
if keepref && length(ref_chn)>1
    signal.data(end+1,:) = ref_signal;
    signal.nbchan = size(signal.data,1);
    signal.chanlocs(end+1).labels = 'CAR';
    signal.chanlocs(end).type = 'REF';
end

exp_endfun('append_online',{'huber_cut',huber_cut});
