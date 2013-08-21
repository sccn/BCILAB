% Local Regression and Likelihood, Figure 10.1.
%
% uses 'module','allcf' rather than 'deriv',[1 1] to get the
% second derivative coefficient
%
% Author: Catherine Loader

load geyser;
a = [0.5 0.75 1];

f=['a','b','c','d'];
for (i=1:3)
  figure('Name', sprintf( 'figure10_1%c', f(i) ) );

  fit = locfit(geyser,'deg',2,'link','ident','kern','gauss','ev','grid','mg',501,'ll',1,'ur',6,'module','allcf','h',a(i));
  x = fit.fit_points.evaluation_points;
  y = lfknots(fit);
  y = y(:,3);
  plot(x,y.^2);
  itgl = [0.5 ones(1,499) 0.5] * (y.^2) / 100
end;

