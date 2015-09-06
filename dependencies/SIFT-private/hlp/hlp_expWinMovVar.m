function [newVar newMean] = hlp_expWinMovVar(varargin)
% this function implements an exponential window moving variance to
% adapt a value based on incoming data
%
% algorithm:
% diff := x - mean 
% incr := alpha * diff 
% mean := mean + incr 
% variance := (1 - alpha) * (variance + diff * incr)

g=arg_define([0 4],varargin, ...
    arg_norep({'values'},mandatory,[],'data values'), ...
    arg_norep({'lastVar','LastVar'},mandatory,[],'previous estimate of the variance'), ...
    arg_norep({'lastMean','LastMean'},mandatory,[],'previous estimate of the mean'), ...
    arg_norep({'numberOfRunsSoFar'},mandatory,[],'number of previous calls to this function'), ...
    arg_norep({'threshold','AbsThresh'},[],[],'Absolute Threshold. Ignore update if absval of the new value exceeds the previous by this factor'), ...
    arg({'mode','Mode'},'cov',{'cov','var'},'Estimate variance or covariance matrix'), ...
    arg({'adaptationHL'},10,[0 Inf],'half-life of exponential win moving average. In frames'), ...
    arg({'updateInterval'},1,[],'num frames between updates'), ...
    arg({'bufferTime'},10,[],'num frames to wait before adapting') ...
    );


% memory factor for exponential window
MEMFACTOR = 2/((g.adaptationHL * 2.8854)+1);
newVar  = g.lastVar;
newMean = g.lastMean;

npts = size(g.values,2);
covmode = strcmpi(g.mode,'cov');
if (g.numberOfRunsSoFar < g.bufferTime || mod(g.numberOfRunsSoFar,g.updateInterval) == 0) && g.numberOfRunsSoFar < Inf
    
    if g.numberOfRunsSoFar < g.bufferTime
        if covmode 
            newVar  = cov(g.values');
        else
            newVar  = var(g.values,[],2);
        end
        newMean = mean(g.values,2);
    else
        % apply exponential-window moving variance

        % get the deviation from the mean
        diff = mean(g.values,2)-g.lastMean;
        % get the increment
        incr = MEMFACTOR * diff;
        % adjust the mean
        newMean = g.lastMean+incr;
        if covmode
            % compute covariance
            newVar  = (1-MEMFACTOR) * (g.lastVar + (diff*incr')/npts);
        else
            % compute variance
            newVar  = (1-MEMFACTOR) * (g.lastVar + (diff.*incr)/npts);
        end
    end
end