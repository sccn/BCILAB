% Local Regression and Likelihood, Figure 5.3.
% Author: Catherine Loader
%
% Density Estimation using poisson regression.
% Stamp thickness data.

load stamp;
n = sum(count);
w = 0.001*n*ones(76,1);
fit = locfit(thick,count,'weights',w,'family','poisson','alpha',[0 0.004]);
figure('Name','fig5_3: Poisson regression density estimation' );
lfplot(fit);
xlabel('Thickness (m.m.)');
ylabel('Density');
