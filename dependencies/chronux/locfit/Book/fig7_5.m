% Local Regression and Likelihood, Figure 7.5.
%
% Author: Catherine Loader

load border;
fit0 = locfit(day,runs,'cens',no,'family','geom','nn',0.7);
fit1 = locfit(day,runs,'weights',ones(265,1)*0.8,'cens',no,'family','geom','nn',0.7);
xev = lfmarg(fit);
xev = xev{1};
y0 = predict(fit0,xev);
y1 = predict(fit1,xev);
figure('Name','fig7_5');
plot(xev,0.8*exp(y1)-exp(y0));
xlabel('Date');
ylabel('Runs');
