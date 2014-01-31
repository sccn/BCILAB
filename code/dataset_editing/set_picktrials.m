function signal = set_picktrials(varargin)
% Select trials whose target variable satisfies a certain condition.
% Result-Set = set_picktrials(Signal, Arguments...)
%
% This function can be used to obtain a subset of the data set, reduced to those epochs that have
% target values matching the supplied criteria.
%
% In:
%   Signal      : Epoched data set with a target variable.
%
%   Arguments...: selection criteria (name-value pairs):
%                    'value': target value or cell array of alternative target values 
%                             -> select by value 
%                    'range': 2-element cell array of the lower and higher end of the range
%                             -> select by range
%                    'rank' : value rank or cell array of alternative values ranks
%                             -> select by rank
%                    currently, only one criterion is supported at a time
%
% Out:
%   Signal  :  Epoched data set which contains only the selected trials.
%
% Notes:
%   This function (in rank-selection mode) is used by methods that need to extract the trials of the
%   first class, second class, etc. from the data.
%
% Examples:
%   % for an epoched data set with target values per epoch, retain only those epochs that
%   % have target value 5 or 7
%   eeg = set_picktrials('Signal',eeg,'ValueSelection',[5 7])
%
%   % for an epoched data set with target values per epoch, retain only those epochs that
%   % have target values between -100 and -50, or between 30 and 40
%   eeg = set_picktrials('Signal',eeg,'RangeSelection',[-100 -50; 30 40])
%  
%   % for an epoched data set with target values per epoch, retain only those epochs that
%   % have the smallest and 2nd smallest target value
%   eeg = set_picktrials('Signal',eeg,'RankSelection',[1 2])
%
% See also:
%   set_picktrials, set_makepos, set_targetmarkers
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01

% set_picktrials_version<1.0> -- for the cache

% sort name-value pairs according to name
if ~exp_beginfun('editing') return; end

declare_properties('independent_channels',true,'independent_trials',false);

arg_define([0 1],varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'valuesel','ValueSelection','value'},[],[],'Set of selected values. Either a single value or a cell array of values.','type','expression','shape','row'), ...
    arg({'rangesel','RangeSelection','range'},[],[],'Range of target values. Two-element cell array of the lower and higher end of the range (inclusive).','type','expression','shape','row'), ...
    arg({'ranksel','RankSelection','rank'},[],[],'Set of selected ranks. Either a single rank value or a cell array of rank values.','type','expression','shape','row'));

if ~isfield(signal,'epoch')
    error('The data set should contain epochs.'); end
if ~isfield(signal.epoch,'target')
    error('The data set does not contain an epoch field named target (denoting the target variable).'); end
if ~isempty(ranksel) + ~isempty(valuesel) + ~isempty(rangesel) ~= 1
    error('Exactly a single criterion is allowed at a time.'); end

variable = vertcat(signal.epoch.target);

if ~isempty(ranksel)
    % do rank-based selection
    if size(variable,2) > 1
        error('Rank-based selection is only supported for one-dimensional target values.'); end
    % rank-transform the variable
    [a,b,variable] = unique(variable); %#ok<ASGLU>
    valuesel = ranksel;
end

if ~isempty(valuesel)
    % do value-based selection
    if ~iscell(valuesel)
        valuesel = {valuesel}; end
    inds = [];
    for i=1:length(valuesel)
        inds = [inds; find(prod(single(variable == repmat(valuesel{i},size(variable,1),1)),2))]; end        
end

if ~isempty(rangesel)
    % do range-based selection
    if ~iscell(rangesel)
        error('Range-based testing requires pairs (two-element cell arrays) of values.'); end
    if ~iscell(rangesel{1})
        rangesel = {rangesel}; end
    inds = [];
    for i=1:length(rangesel)
        if length(rangesel{i}) ~= 2
            error('Range-based testing requires pairs (two-element cell arrays) of values.'); end
        inds = [inds; intersect(find(prod(single(variable >= repmat(rangesel{i}{1},size(variable,1),1)),2)), ...
                                find(prod(single(variable <= repmat(rangesel{i}{2},size(variable,1),1)),2)))]; 
    end
end

if isempty(inds)
    error('BCILAB:set_picktrial:no_trials','This data set contains no trials for one of your target classes: please check whether your target markers are correct.\n\nAlso note that while the entire recording might have trials for each class, the subset used for training may not, for example when doing a block-wise cross-validation on a block experiment design that has no repetitions.');
else
    signal = exp_eval(set_selepos(signal,inds));
end

exp_endfun;
