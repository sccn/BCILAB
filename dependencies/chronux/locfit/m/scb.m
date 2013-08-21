function z=scb(x,y,varargin)

% Simultaneous Confidence Bands
%
% Example:
%   load ethanol;
%   z = scb(E,NOx,'h',0.5);
%
% result (z) is a matrix with four columns: evaluation points,
% fitted values, lower confidence limit, upper confidence limit.
% Most locfit arguments should work.

fit = locfit(x,y,'ev','grid','mg',20,varargin{:});
kap = kappa0(x,y,varargin{:});
cb = predict(fit,'fitp','band','g','kappa',kap);
z = [fit.fit_points.evaluation_points' cb{1} cb{3}];

return;
