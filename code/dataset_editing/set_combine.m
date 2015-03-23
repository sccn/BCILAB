function result = set_combine(varargin)
% Combine data sets into a new multi-class data set (with multiple target values).
% Result = set_combine(Set-#1, Set-#2, ...)
%
% Combines the specified epoched data sets into a new data set, which has a target value per epoch,
% determined by the index of the data set from which the respective epoch was taken (index as in
% set_combine's argument list).
%
% When data sets are recorded or imported, they are usually not annotated with a "target variable"
% per trial, which is however required to calibrate predictive models so that they can make
% predictions (of that target variable). The target variable is, thus, a variable
% (categorical/real/multivariate) that can be assigned to a trial, and paradigms learn to predict
% this variable for a given trial. Usually, this variable is derived from events that have been
% recorded together with the data set, and the most direct approach is to map certain event types to
% certain target values. This can be accomplished via the 'events' option in most paradigms and
% filters/flt_pipeline, or directly via dataset_ops/set_makepos.
%
% If the mapping from event types to target variable is very complex, a great deal of event
% rewriting may be necessary until one of these methods can be applied, and the alternative is to
% handle each condition in a separate data set, using EEGLAB functions (deleting all trials that are
% not in the condition, depending on the events in the trial), and then doing a final merge of all
% those data sets using set_combine, where the target variable for all trials that came from the
% first specified set is taken to be 1, for all trials that came from the second set is set to 2,
% etc.
%
% In:
%   Set-#k  : epoched data set #k
%
% Out:
%   Result  : the combined data set, containing all trials of the given sets, but with the target 
%             variable of trials from Set-#i assigned i.
%
% Notes:
%   The epochs of both classes must have the same number of time points. Parameters cannot be passed
%   by name to this function.
%
% Examples:
%   % concatenate the epochs of the dataset condition_A and the dataset condition_B, and assign
%   % the field eeg.epoch.target such that all epochs of the first set receive target value 1 and all
%   % epochs of the second set receive target value 2.
%   eeg = set_combine(condition_A,condition_B)
%
% See also:
%   set_merge
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-31
dp;

% set_combine_version<1.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('name','CombineSets','independent_channels',true,'independent_trials',true);

if isempty(varargin)
    result = []; 
else
    % do a regular merge
    [result,indexmap] = exp_eval(set_merge(varargin{:}));
    if ~isempty(result.epoch)
        % but set each epoch's target value to the k of the respective combined set
        indexmap = num2cell(indexmap);
        if length(result.epoch) < length(indexmap)
            result.epoch(length(indexmap)).target = []; end
        [result.epoch.target] = indexmap{:};
    else
        disp_once('Warning: the dataset after set_merge() does not appear to be epoched.'); 
    end
end

% take over the online expression of the merge result
exp_endfun;
