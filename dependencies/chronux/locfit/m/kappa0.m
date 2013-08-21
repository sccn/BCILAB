function kap=kappa0(x,y,varargin)

% Compute the constants for `tube-formula' based simultaneous
% confidence bands.
%
% Works for regression models only. Density estimation problems
% should be converted to counts, and use poisson regression
% 'family','poisson'.
%
% Essentially, this is a front-end to locfit, and so all optional
% arguments to locfit (eg, smoothing parameters) can be provided.
%
% To compute (or plot) the confidence bands, provide the output
% of the kappa0() function as the 'kappa' argument to a
% predict() or lfband() call.
%
%
% Example:
%
% load ethanol;
% fit = locfit(E,NOx,'alpha',0.5)
% kap = kappa0(E,NOx,'alpha',0.5)  % give same arguments!
% lfplot(fit)
% lfband(fit,'kappa',kap)     % plot the simultaneous bands
% z = predict(fit,[0.6 0.7 0.8]','kappa',kap,'band','g')
% z{3}                        % evaluate the bands.

fit = locfit(x,y,'module','kappa','ev','grid','mg',20,varargin{:});
z = fit.fit_points.kappa;
d = size(fit.data.x,2);
kap = z(1:(d+1));

return;
