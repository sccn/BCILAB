function [limits state] = hlp_adaptMinMaxLimits(varargin)
% adapt min and max limits using an exponential window moving average. The
% moving average is calculated in place based on new data values and
% previous min/max limits.

g = arg_define(varargin, ...
        arg_norep({'values'},[],[],'data values. can be scalar, vector or matrix'), ...
        arg_nogui({'state'},[],[],'current state of the adaptation'), ...
        arg_sub({'adaptOpts'},{},@hlp_scaleLimits,'Adaptation options'), ...
        arg({'reset'},false,[],'Reset adaptation state') ...
        );

state = g.state;

if g.reset || isempty(state)
    % initialize state
    state.lastMin      = min(g.values(:));
    state.lastMax      = max(g.values(:));
    state.numRunsSoFar = 0;
end

if isfield(g.adaptOpts,'arg_direct')
    g.adaptOpts = rmfield(g.adaptOpts,'arg_direct');
end

% adapt limits using exponential window moving average
[state.lastMin state.lastMax] = hlp_scaleLimits(g.values,state.lastMin, ...
                                state.lastMax,state.numRunsSoFar,g.adaptOpts);
state.numRunsSoFar = state.numRunsSoFar + 1;

limits = [state.lastMin, state.lastMax];