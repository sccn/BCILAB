% Local Regression and Likelihood, Figure 6.1.
%
% Derivative (local slope) estimation, for the Old Faithful Geyser Data.
% The 'deriv' argument specifies derivative estimation,
%   'deriv',1        First-order derivative.
%   'deriv',[1 1]    Second-order derivative.
%   'deriv',2        For bivariate fits, partial deriv. wrt second variable.
%   'deriv',[1 2]    Mixed second-order derivative.
%
% Density estimation is done on the log-scale. That is, the estimate
% is of g(x) = log(f(x)), where f(x) is the density.
%
% The relation between derivatives is therefore
%     f'(x) = f(x)g'(x) = g'(x)exp(g(x)).
% To estimate f'(x), we must estimate g(x) and g'(x) (fit1 and fit2 below),
% evaluate on a grid of points (p1 and p2), and apply the back-transformation.
%
% Disclaimer: I don't consider derivative estimation from noisy data
% to be a well-defined problem. Use at your own risk.
%
% Author: Catherine Loader
%
% NEED: m argument passed to lfmarg().

load geyser;
fit1 = locfit(geyser,'alpha',[0.1 0.6],'ll',1,'ur',6);
fit2 = locfit(geyser,'alpha',[0.1 0.6],'ll',1,'ur',6,'deriv',1);
z = lfmarg(fit1);
p1 = predict(fit1,z);
p2 = predict(fit2,z);
figure('Name','fig6_1: slope estimation: Old faithful data' );
plot(z{1},p2.*exp(p1));
xlabel('Eruption Duration (Minutes)');
ylabel('Density Derivative');
