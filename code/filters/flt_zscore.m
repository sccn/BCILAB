function [signal state] = flt_zscore(varargin)
% Exponential window moving zscore

if ~exp_beginfun('filter'), return; end

declare_properties('name', 'AdaptiveZscore', 'experimental',true, 'experimental',true, 'cannot_follow','set_makepos', 'independent_trials',true, 'independent_channels',true);

state = [];
arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg_sub({'adaptOpts'},{'adaptationHL',100},@hlp_expWinMovVar,'Adaptation options'), ...
    arg({'doreset','reset'},false,[],'Reset adaptation state'), ...
    arg({'verb','Verbosity'},false,[],'Verbose output'), ...
    arg_norep({'state','State'},unassigned));

if doreset || isempty(state)
    % initialize state
    state.lastVar      = 1;
    state.lastMean     = 0;
    state.numRunsSoFar = 0;
end

% adaptive z-scoring
% adapt limits using exponential window moving average
[state.lastVar state.lastMean] = hlp_expWinMovVar('values',signal.data,  ...
                                                  'lastVar',state.lastVar,      ...
                                                  'lastMean',state.lastMean,     ...
                                                  'numberOfRunsSoFar',state.numRunsSoFar, ...
                                                  rmfield(adaptOpts,'arg_direct'));
state.numRunsSoFar = state.numRunsSoFar + 1;

signal.data = (signal.data - state.lastMean)./ sqrt(state.lastVar);

if verb
    fprintf('data:%s\nmean:%s\nvar:%s\n\n',hlp_tostring(signal.data),hlp_tostring(state.lastMean),hlp_tostring(state.lastVar));
end

exp_endfun;
