% Local Regression and Likelihood, Figure 4.3.
% Author: Catherine Loader
%
% Residual plots for local likelihood.
% Henderson-Shepherd Mortality Data.

load morths;
fit = locfit(age,deaths,'weights',n,'family','binomial','alpha',0.5);

figure('Name','fig4_3a: Residual plots for local likelihood');
plot(age,residuals(fit,'dev'),'.-');
xlabel('Age');
ylabel('Residual');
title('Deviance Residuals');
hold on;
plot([min(age) max(age)],[0 0],':');
hold off;

figure('Name','fig4_3b: Residual plots for local likelihood');
plot(age,residuals(fit,'pear'),'.-');
xlabel('Age');
ylabel('Residual');
title('Pearson Residuals');
hold on;
plot([min(age) max(age)],[0 0],':');
hold off;

figure('Name','fig4_3c: Residual plots for local likelihood');

plot(age,residuals(fit,'raw'),'.-');
xlabel('Age');
ylabel('Residual');
title('Raw (response) Residuals');
hold on;
plot([min(age) max(age)],[0 0],':');
hold off;

figure('Name','fig4_3d: Residual plots for local likelihood');
plot(age,residuals(fit,'ldot'),'.-');
xlabel('Age');
ylabel('Residual');
title('ldot Residuals');
hold on;
plot([min(age) max(age)],[0 0],':');
hold off;

