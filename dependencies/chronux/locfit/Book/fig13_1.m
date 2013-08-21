% Local Regression and Likelihood, Figure 13.1
%
% Author: Catherine Loader

load mmsamp;

figure('Name','fig13_1a');
plot(x,y,'o');
hold on;

fit = locfit(x,y,'deg',1,'kern','minmax','alpha',4000,'ev','grid','mg',200,'ll',0,'ur',1);
xev = fit.fit_points.evaluation_points';
yev = predict(fit,'fitp');
plot(xev,yev);

fit = locfit(x,y,'deg',1,'h',0.05,'ev','grid','mg',100,'ll',0,'ur',1);
xev = fit.fit_points.evaluation_points';
yev = predict(fit,'fitp');
plot(xev,yev,'red');

xx = 0:0.01:1;
yy = 2-5*xx+5*exp(-(20*xx-10).*(20*xx-10));
plot(xx,yy,'green');

legend('Data','Minimax','Constant h','True Mean');

hold off;
