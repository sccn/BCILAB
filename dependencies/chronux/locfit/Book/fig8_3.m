% Local Regression and Likelihood, Figure 8.3.
%
% Discrimination/Classification, iris data for different
% smooting paramters.
%
% Note that the iris.mat file contains the full iris dataset;
% only Versicolor and Virginica are used in this example.
%
% The `Species' variable contains species names. `Specn' has
% them numerically, Setosa=1, Versicolor=2, Virginica=3.
%
% Author: Catherine Loader
%
% Need: get contour plot correct.
% Need: distinguish colors in plot.

load iris;
a = (2:9)/10;
z = zeros(size(a));

u = find(Specn >= 2);
pw = PetalWid(u);
pl = PetalLen(u);
y = (Specn(u)==3);

for i = 1:length(a)
  fit = locfit([pw pl],y,'deg',1,'alpha',a(i),'ev','cros','scale',0,'family','binomial');
  fv = fitted(fit);
  tb = tabulate(10*y+(fv>=0.5))
  z(i) = sum(y == (fv<=0.5));
end;

[a; z]

fit = locfit([pw pl],y,'deg',1,'scale',0,'family','binomial');
figure('Name','fig8_3: iris classification');
lfplot(fit,'contour');

