% Local Regression and Likelihood, Figure 6.4.
% Author: Catherine Loader
%
% Local smooth of CO2 dataset. Estimate the main trend,
% then use periodic smoothing of the residuals to estimate
% the annual effect. Add main trend and periodic components
% to get overall smooth.
%
% A periodic smooth is specified by 'style','a'.
% Note that year+month/12 scales the predictor to have a period
% of 1. The 'scale' argument to locfit() is period/(2*pi).

load co2;
fit1 = locfit(year+month/12,co2,'alpha',0.5,'deg',1);
res = residuals(fit1);
fit2 = locfit(year+month/12,res,'alpha',[0 2],'style','a','scale',1/(2*pi));
f1 = fitted(fit1);
f2 = fitted(fit2);
figure('Name','fig6_4: CO2 dataset local smoothing' );
plot(year+month/12,f1+f2);
