% Local Regression and Likelihood, Figure 7.6.
%
% Estimating the mean survival time using a local weibull
% (transformed Gamma) model.

load heart;
z = 0.625;
ts = (surv+0.5).^z;
fit = locfit(age,ts,'cens',cens,'family','gamma','nn',0.8);
i = find(cens==0);
figure('Name','fig7_6: mean survival time local weibull');
plot(age(i),log(surv(i)+0.5),'o');
hold on;
i = find(cens==1);
plot(age(i),log(surv(i)+0.5),'r+');
xlabel('Age at Transplant (Years)');
ylabel('log(0.5+Survival Time (Days))');

xev = lfmarg(fit);
y = predict(fit,xev);
plot(xev{1},log(exp(y./z)*gamma(1+1./z)+0.5));
hold off;

