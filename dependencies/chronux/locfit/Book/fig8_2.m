% Local Regression and Likelihood, Figure 8.2.
%
% Discrimination/Classification, simple example using
% density estimation.
%
% First, compute density estimates fit0, fit1 ('family','rate'
% - output is in events per unit area) for each class in the
% training sample. The ratio fit1/(fit1+fit0) estimates the
% posterior probability that an observation comes from population 1.
%
% plotting the classification boundary is slightly tricky - it depends
% on both fits, so lfplot() can't be used. Instead, both fits must be
% evaluated on the same grid of values, which is then used to make a
% contour plot.
%
% Author: Catherine Loader

load cltrain;
u0 = find(y==0);
u1 = find(y==1);
fit0 = locfit([x1(u0) x2(u0)],y(u0),'family','rate','scale',0);
fit1 = locfit([x1(u1) x2(u1)],y(u1),'family','rate','scale',0);

v0 = -3+6*(0:50)'/50;
v1 = -2.2+4.2*(0:49)'/49;
% predict returns log(rate)
z = predict(fit0,{v0 v1})-predict(fit1,{v0 v1});
z = reshape(z,51,50);
figure('Name','fig8_2: classification');
contour(v0,v1,z',[0 0]);
hold on;
plot(x1(u0),x2(u0),'.');
plot(x1(u1),x2(u1),'.','color','red');
hold off;

p0 = predict(fit0,[x1 x2]);
p1 = predict(fit1,[x1 x2]);
py = (p1 > p0);
disp('Classification table for training data');
tabulate(10*y+py);

load cltest;
p0 = predict(fit0,[x1 x2]);
p1 = predict(fit1,[x1 x2]);
py = (p1 > p0);
disp('Classification table for test data');
tabulate(10*y+py);

