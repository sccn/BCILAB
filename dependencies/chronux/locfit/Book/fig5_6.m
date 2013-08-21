% Local Regression and Likelihood, Figure 5.6.
% Author: Catherine Loader
%
% Cumulative Distribution Function
%
% This function uses cdfplot() from Matlab's add-on statistics toolbox.
%
% NEED:
%   predict should backtr by default?

load geyser;
fit = locfit(geyser,'nn',0.1,'h',1.2,'renorm',1);
x = (1:0.01:6)';
z = predict(fit,x);
figure('Name','fig5_6: Cumulative distribution function' );
plot(x,0.01*cumsum(exp(z)));
xlabel('Eruption Duration (Minutes)');
ylabel('Cumulative Distribution Function');
hold on;
cdfplot(geyser);
hold off;
