function signal = flt_selchans(varargin)
% Selects a subset of channels from the given data set.
% Signal = flt_selchans(Signal, Channels, OrderPreservation, RemoveSelection, FindClosest)
%
% Channel (or sensor) selection is a simple and effective method to constrain the complexity (and
% thus shorten computation time and/or improve robustness) of later stages in a paradigm. Sometimes,
% it is also employed to approximately restrict a paradigm to a sub-portion of brain processes (by
% using only channels directly above the brain region of interest), but this is not guaranteed to
% have the desired effect, since eletromagnetic signals emitted from any point in the brain are
% practically captured by every sensor (due to volume conduction). Other uses of channel selection
% are to exclude bad channels in a faulty recording or to simulate the behavior of a paradigm
% running on a subset of the sensors (e.g., for cost reduction purposes).
% 
% In:
%   Signal    : Data set
%
%   Channels  : channel indices or names to select
%
%   OrderPreservation : Output channel order. The result will have its channels either in the order of the 
%                       input set (if 'dataset-order') or in the order of the query list (if 'query-order').
%                       (default: 'query-order')
%
%   RemoveSelection : Remove selected channels. (default: false)
%
%   FindClosest : Find closest channels. This is for cases where the requested channels are not in
%                 the set. (default: false)
%
% Out:
%   Signal    : The original data set restricted to the selected channels (as far as they are 
%               contained)
%
% Examples:
%   % select only the channels C3, C4 and Cz
%   eeg = flt_selchans(eeg,{'C3','C4','Cz'})
%
%   % select channels 1:32
%   eeg = flt_selchans(eeg,1:32)
%
%   % reverse channel order
%   eeg = flt_selchans(eeg,eeg.nbchan:-1:1)
%
%   % retain all channels (i.e., do nothing)
%   eeg = flt_selchans(eeg,[])
%
%   % select a group of channels but keep them in the order in which they were in the original data set
%   eeg = flt_selchans(eeg,{'AFz','Fz','Fpz','F1'},'dataset-order')
%
% See also:
%   flt_seltypes
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-17

% flt_selchans_version<1.11> -- for the cache

if ~exp_beginfun('filter') return; end

% used as a tool to select channel subsets before these ops are applied
declare_properties('name',{'ChannelSelection','channels'}, 'precedes',{'flt_laplace','flt_ica','flt_reref'}, 'independent_trials',true, 'independent_channels',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'channels','Channels'}, [], [], 'Cell array of channel names to retain.','type','cellstr','shape','row'), ...
    arg({'orderPreservation','OrderPreservation'}, 'query-order', {'query-order','dataset-order'}, 'Output channel order. The result will have its channels either in the order of the input set or in the order of the query list.'), ...
    arg({'remove_selection','RemoveSelection'},false,[],'Remove selected channels.'), ...
    arg({'find_closest','FindClosest'},false,[],'Find closest channels. This is for cases where the requested channels are not in the set.'));

% determine channel indices to retain
if find_closest
    if ~iscellstr(channels)
        error('For distance-based channel matching the given channels should be a cell array of labels.'); end
    tmplocs = hlp_microcache('matchchan',@set_infer_chanlocs,channels);
    subset = hlp_microcache('matchchan',@eeg_matchchans,signal.chanlocs,tmplocs,'noplot');
else
    subset = set_chanid(signal,channels);
end

% optionally invert selection
if remove_selection
    tmp = true(1,signal.nbchan);
    tmp(subset) = false;
    subset = find(tmp);
end

% handle order preservation
if strcmp(orderPreservation,'dataset-order')
    subset = sort(subset);
elseif ~strcmp(orderPreservation,'query-order')
    error(['Unknown order requested: ' orderPreservation]);
end

if ~isequal(subset,1:signal.nbchan)
    % update .data
    signal.data = signal.data(subset,:,:,:,:,:,:,:);
    % update .chanlocs and .nbchan
    signal.chanlocs = signal.chanlocs(subset);
    signal.nbchan = size(signal.data,1);
    % reset any ICA parameters (should be recomputed)
    signal.icachansind = [];
    signal.icaweights = [];
    signal.icasphere = [];
    signal.icawinv = [];
end

exp_endfun;
