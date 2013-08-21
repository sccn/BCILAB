% Local Regression and Likelihood, Figure 5.4.
% Author: Catherine Loader
%
% Bivariate Density Estimation.

load trimod;
fit = locfit([x0 x1],'nn',0.35);
figure('Name','fig5_4: Bivariate density estimation' );
lfplot(fit);
