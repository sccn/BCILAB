% Local Regression and Likelihood, Figure 3.1.
%
% Bivariate Local Regression.

load ethanol;

fit = locfit([E C],NOx,'alpha',0.5);
figure('Name','fig3_1: Bivariate local Regression');
lfplot(fit);
