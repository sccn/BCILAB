% Local Regression and Likelihood, Figure 8.1.
%
% Discrimination/Classification, simple example using
% logistic regression.
%
% Author: Catherine Loader
%
% NEED: level=0.5 for contour plot.

load cltrain;
fit = locfit([x1 x2],y,'family','binomial','scale',0);
figure('Name','fig8_1: logistic regression');
lfplot(fit,'contour');
hold on;
u = find(y==0);
plot(x1(u),x2(u),'.');
u = find(y==1);
plot(x1(u),x2(u),'.','color','red');
hold off;

py = fitted(fit)>0.5;
disp('Classification table for training data');
tabulate(10*y+py);

load cltest;
py = predict(fit,[x1 x2])>0;
disp('Classification table for test data');
tabulate(10*y+py);
