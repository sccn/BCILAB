% Local Regression and Likelihood, Figure 11.5.
% Author: Catherine Loader
%
% Local Adaptive Smooth - Dopler example.
% Adaptive fitting using Katkovnik (et al) ICI criterion.
%
% NEED: Improve plot - curve is hard to see.

x = (0:2047)'/2047;
m = 20*sqrt(x.*(1-x)) .* sin(2*pi*1.05 ./ (x+0.05));
y = m + normrnd(0,1,2048,1);
fit = locfit(x,y,'pen',1.1,'acri','ici','maxk',500);
figure('Name','fig11_5: local adaptive smooth: dopler example');
lfplot(fit);
