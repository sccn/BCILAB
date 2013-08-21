% Local Regression and Likelihood, Figure 2.2.
%
% Using different smoothing parameters.

load ethanol;
fit = locfit(E,NOx,'alpha',0.5);
figure('Name', 'figure2_2: Smoothing parameters' );
lfplot(fit);
