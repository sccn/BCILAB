% Local Regression and Likelihood, Figure 6.6.
% Author: Catherine Loader
%
% Penny data, estimating jump size using left and
% right smooths.
%

load penny;
midp = (1945:1988)'+0.5;
fitl = locfit(year,thickness,'style','l','alpha',[0 10],'ev',midp,'deg',1);
fitr = locfit(year,thickness,'style','r','alpha',[0 10],'ev',midp,'deg',1);
pl = predict(fitl);
pr = predict(fitr);
figure('Name','fig6_6: Estimating jump size' );
plot(midp,(pr-pl).*(pr-pl),'.-');
xlabel('Year');
ylabel('\Delta^2(Year)');
