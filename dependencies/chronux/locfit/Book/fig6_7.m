% Local Regression and Likelihood, Figure 6.6.
% Author: Catherine Loader
%
% Penny data, piecewise smoothing.
%

load penny;
u = find(year <= 1958);
fit1 = locfit(year(u),thickness(u),'alpha',[0 10],'deg',1);
u = find((year>1958) & (year <= 1974));
fit2 = locfit(year(u),thickness(u),'alpha',[0 10],'deg',1);
u = find(year > 1974);
fit3 = locfit(year(u),thickness(u),'alpha',[0 10],'deg',1);
figure('Name','fig6_7: Piecewise smoothing of penny data' );
lfplot(fit1);
hold on;
lfplot(fit2);
hold on;
lfplot(fit3);
