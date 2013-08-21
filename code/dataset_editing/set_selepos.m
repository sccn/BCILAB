function signal = set_selepos(varargin)
% Select a subset of epochs from a data set.
% Signal = set_selepos(Signal, EpochIndices)
%
% In:
%   Signal       : epoched EEGLAB data set from which epochs shall be selected
%
%   EpochIndices : indices of the epochs that should be retained
% 
% Out:
%   Signal       : newly created data set that contains only the selected epochs.
%
% Examples:
%   % for an epoched data set, select epochs 50-100
%   eeg = set_selepos(eeg,50:100)
%
% See also:
%   set_selinterval, set_partition, set_makepos
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01

% set_selepos_version<1.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('independent_channels',true,'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'epos','EpochIndices'},[],[],'Indices of retained epochs.','shape','row'));

signal = pop_select(signal,'trial',epos,'sorttrial','off');

exp_endfun;
