function [newVar newMean state] = hlp_adaptExpWinMovVar(varargin)
% estimate exponential window moving variance (and mean). The
% moving variance is calculated in place based on new data values and
% previous data.

g = arg_define(varargin, ...
        arg_norep({'values'},[],[],'data values. can be scalar, vector or matrix'), ...
        arg_nogui({'instate','State'},[],[],'state'), ...
        arg_sub({'adaptOpts'},{},@hlp_scaleLimits,'Adaptation options'), ...
        arg({'reset'},false,[],'Reset adaptation state') ...
        );

state = g.instate;

if g.reset || ~isfield(g,'instate') || isempty(state)
    % initialize state
    state.lastVar      = 1;
    state.lastMean     = 0;
    state.numRunsSoFar = 0;
end

% adapt limits using exponential window moving average
[state.lastVar state.lastMean] = hlp_expWinMovVar(g.values,state.lastVar, state.lastMean, ...
                                     state.numRunsSoFar,rmfield(g.adaptOpts,'arg_direct'));
state.numRunsSoFar = state.numRunsSoFar + 1;

newVar  = state.lastVar;
newMean = state.lastMean;