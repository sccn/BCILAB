% Local Regression and Likelihood, Figure 5.5.
% Author: Catherine Loader
%
% Bivariate Density Estimation.

load trimod;
fit = locfit([x0 x1],'nn',0.35);
emp = sort(fitted(fit));
v = emp(floor([0.05 0.5]*225));
figure('Name','fig5_5: Bivariate density estimation' );
lfplot(fit,'contour',v);
