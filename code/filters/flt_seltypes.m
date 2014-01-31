function signal = flt_seltypes(varargin)
% Selects a subset of channels based on their type from the given data set.
% Signal = flt_seltypes(Signal, Types)
%
% This function allows to select a subset of channels in a given data set that matches a particular
% channel type. The channel type is stored in the field .chanlocs(i).type (for channel i), and is 
% also auto-inferred by io_loadset roughly based on typical channel names.
% 
% In:
%   Signal    : Data set
%
%   ChannelTypes     : channel type string or cell array of channel types to select
%
% Out:
%   Signal    : The original data set restricted to the selected channels (as far as they are 
%               contained)
%
% Examples:
%   % select only the channels C3, C4 and Cz
%   eeg = flt_seltypes(data,'EEG')
%
%   % select channels 1:32
%   eye = flt_seltypes(eeg,'EOG')
%
% See also:
%   flt_selchans, io_loadset
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-11-23

% flt_seltypes_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

% used as a tool to select channel subsets before these ops are applied
declare_properties('name',{'TypeSelection','types'}, 'precedes',{'flt_laplace','flt_ica','flt_reref'}, 'independent_trials',true, 'independent_channels',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'chantypes','ChannelTypes'}, [], [], 'Cell array of channel types to retain.','type','cellstr','shape','row'));

if ischar(chantypes)
    chantypes = {chantypes}; end

retain_indices = [];
for t = chantypes
    retain_indices = [retain_indices find(strcmp({signal.chanlocs.type},t{1}))]; end %#ok<AGROW>
signal.data = signal.data(sort(retain_indices),:,:,:,:,:,:,:);
signal.chanlocs = signal.chanlocs(sort(retain_indices));
signal.nbchan = size(signal.data,1);

exp_endfun;
