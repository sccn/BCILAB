function signal = flt_reref(varargin)
% Re-references the data to a new (set of) channel(s) or the average of all channels.
% Signal = flt_reref(Signal, ReferenceChannels, ExcludeChannels, KeepReference)
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
% See also:
%   pop_reref
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% flt_reref_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

declare_properties('name',{'Rereferencing','ref'}, 'independent_channels',false,'independent_trials',true);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'chn','ReferenceChannels'}, [], [], 'Cell array of reference channels. The signal data is be referenced to these, defaults to average reference if empty.','type','cellstr','shape','row'),...
    arg({'exclude','ExcludeChannels'}, [], [], 'Cell array of channels to exclude.','type','cellstr','shape','row'),...
    arg({'keepref','KeepReference'}, false, [], 'Keep the reference channel.'));

signal = pop_reref(signal,set_chanid(signal,chn),'exclude',fastif(isempty(exclude),[],set_chanid(signal,exclude)),'keepref',fastif(keepref,'on','off')); %#ok<*NODEF>

exp_endfun;
