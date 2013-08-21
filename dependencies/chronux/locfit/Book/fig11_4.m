% Local Regression and Likelihood, Figure 11.4.
% Author: Catherine Loader
%
% Local Adaptive Smooth - Dopler example.
%
% NEED: Improve plot - curve is hard to see.

x = (0:2047)'/2047;
m = 20*sqrt(x.*(1-x)) .* sin(2*pi*1.05 ./ (x+0.05));
y = m + normrnd(0,1,2048,1);
fit = locfit(x,y,'pen',4,'acri','cp','maxk',500);
figure('Name','fig11_4: local adaptive smooth: dopler example');
lfplot(fit,'red');
