% Local Regression and Likelihood, Figure 4.2.
% Author: Catherine Loader
%
% Local Likelihood (Logistic Regression), for the
% Henderson-Shepherd Mortality Data.

load morths;
fit = locfit(age,deaths,'weights',n,'family','binomial','alpha',0.5);
figure('Name','fig4_2: Logistic Regression for the Henderson-Shepherd Mortality Data.');
lfplot(fit);
lfband(fit);
