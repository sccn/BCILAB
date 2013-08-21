% Local Regression and Likelihood, Figure 7.3.
%
% Censored local regression.
%
% Author: Catherine Loader
%
% NEEDS: Kaplan-Meier based estimate. Identify curves and cens/uncens. data

load heart;
fit = locfit(age,log(0.5+surv));
figure('Name','fig7_3: censored local regression Kaplan-Meier based;' );
lfplot(fit);
xlabel('Age at Transplant (years)');
ylabel('0.5+Survival Time (Days)');

fit = lf_censor(age,log(0.5+surv),cens);
hold on;
lfplot(fit,'nodata');
hold off;
