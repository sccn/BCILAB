function y = logistic(x,transWidth,inflectPoint,max,transThresh)
% evaluate the function y=f(x) where f(x) is the logistic (sigmoid) pdf
%
% transWidth: transition width (samples) -- we define this as the width W
%                       for which
%                       y(inflectPoint+W)   = (1-q)*max{y} and
%                       y(inflectPoint-W)   = q*max{y}
%             in other words, y is within q*100% of its asymptote at this
%             point
% monotonically decreasing. If positive, then increasing)
% inflectPoint: inflection point
% max: maximum value
%
% (C) Tim Mullen, 2011. SCCN/INC UCSD


if nargin<5
    q = 0.1; % percent of max for transition width
else
    q = transThresh;
end

beta = transWidth/log(1/q - 1);

y = max./(1+exp(-(1/beta).*(x-inflectPoint)));