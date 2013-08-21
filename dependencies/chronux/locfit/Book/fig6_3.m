% Local Regression and Likelihood, Figure 6.3.
% Author: Catherine Loader
%
% Local smooth of CO2 dataset. Estimate the main trend,
% then use periodic smoothing of the residuals to estimate
% the annual effect.
%
% A periodic smooth is specified by 'style','a'.
% Note that year+month/12 scales the predictor to have a period
% of 1. The 'scale' argument to locfit() is period/(2*pi).

load co2;
fit1 = locfit(year+month/12,co2,'alpha',0.5,'deg',1);
res = residuals(fit1);
fit2 = locfit(year+month/12,res,'alpha',[0 2],'style','a','scale',1/(2*pi));
figure('Name','fig6_3: CO2 Dataset: Local smoothing' );
lfplot(fit2,'nodata');
