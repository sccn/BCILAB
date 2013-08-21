% Local Regression and Likelihood, Figure 4.1.
% Author: Catherine Loader
%
% Local Likelihood (Poisson Regression).

load mine;
fit = locfit(extrp,frac,'family','poisson','deg',1,'alpha',0.6);
figure('Name','fig4_1: Poisson Regression');
lfplot(fit);
lfband(fit);
