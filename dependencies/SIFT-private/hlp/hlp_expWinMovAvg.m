function newVals = hlp_expWinMovAvg(varargin)

% this function implements an exponential window moving average to
% adapt a value based on incoming data

g=arg_define([0 3],varargin, ...
    arg_norep({'values'},mandatory,[],'data values'), ...
    arg_norep({'lastVal','LastValue'},mandatory,[],'last value'), ...
    arg_norep({'numberOfRunsSoFar'},mandatory,[],'number of previous calls to this function'), ...
    arg_norep({'threshold','AbsThresh'},[],[],'Absolute Threshold. Ignore update if absval of the new value exceeds the previous by this factor'), ...
    arg({'adaptationHL'},10,[0 Inf],'half-life of exponential win moving average. In frames'), ...
    arg({'updateInterval'},1,[],'num frames between updates'), ...
    arg({'bufferTime'},10,[],'num frames to wait before adapting') ...
    );


% memory factor for exponential window moving average
MEMFACTOR = 2/((g.adaptationHL * 2.8854)+1);

newVals = g.lastVal;

if (g.numberOfRunsSoFar < g.bufferTime || mod(g.numberOfRunsSoFar,g.updateInterval) == 0) && g.numberOfRunsSoFar < Inf
    

    if g.numberOfRunsSoFar < g.bufferTime
        newVals = g.values;
    else
        % apply exponential-window moving average to calculate new value
        newVals  = MEMFACTOR * g.values(:) + (1-MEMFACTOR) * g.lastVal;
    end
    
%     if ~isempty(g.threshold) && g.lastVal>0 ...
%         && g.values(:) < g.lastMin*g.threshold(1) || max(g.values(:)) > g.lastMax*g.threshold(2))
%         newMin = g.lastMin;
%         newMax = g.lastMax;
%         return;
%     end
    
end