function yhat=lfsmooth(varargin)
%
% a simple interface to locfit.
% output is a vector of smoothed values, at each data point.
% all locfit options, except evaluation structures, are valid.
%
% Example, to smooth a time series of observations,
%
% t = (1:100)';
% y = 2*sin(t/10) + normrnd(0,1,100,1);
% plot(t,y,'.');
% hold on;
% plot(t,lfsmooth(t,y,'nn',0.5));
% hold off;
%

% Minimal input validation    
if nargin < 1
   error( 'At least one input argument required' );
end

fit = locfit(x,varargin{:},'module','simple');
yhat = fit.fit_points.fitted_values;

return;
