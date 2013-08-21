% Local Regression and Likelihood, Figure 6.5.
% Author: Catherine Loader
%
% Local smooth of CO2 dataset.
% Use bivariate fit to capture the two trends.
%
% bonus plot: fitted values vs year.

load co2;
fit = locfit([month year+month/12],co2,'alpha',0.2,'style','an','scale',[6/pi,10]);

figure('Name','fig6_5a: Local smoothing');
lfplot(fit,'contour');
xlabel('Month');
ylabel('Year');

figure('Name','fig6_5b: Local smoothing');
plot(year+month/12,fitted(fit));
xlabel('Year');
ylabel('CO_2');

