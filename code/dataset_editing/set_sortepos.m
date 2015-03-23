function [signal,inds] = set_sortepos(varargin)
% Sort epochs in a given data set by their latency, for chronological cross-validation.
% [Signal,Indices] = set_sortepos(Signal)
%
% In:
%   Signal : epoched data set for which epochs shall be sorted
%
% Out:
%   Signal : epoched data set with sorted epochs
%
% Note:
%   Sorts according to the EEG.epoch.latency field.
%
% Examples:
%   % for an epoched data set, sort the epochs by time (of their time-locking event)
%   eeg = set_sortepos(eeg)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-31

% set_sortepos_version<1.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('independent_channels',true,'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}));

if isfield(signal,'epoch') %#ok<NODEF>
    if ~isempty(signal.epoch)
        if ~isfield(signal.epoch,'latency')
            error('Your signal is lacking the required .epoch.latency field.'); end
        [x,inds] = sort([signal.epoch.latency]); %#ok<ASGLU>
        signal = exp_eval(set_selepos(signal,inds));
    end
else
    error('The given signal parameter is lacking the required .epoch field. Not a valid signal.');
end

exp_endfun;
