% Local Regression and Likelihood, Figure 4.4.
% Author: Catherine Loader
%
% AIC plot for a local poisson regression.

load mine;
a = (0.4:0.05:1)';
figure('Name','fig4_4: AIC plot for local poisson regression' );
aicplot(a,extrp,frac,'family','poisson','deg',1);
