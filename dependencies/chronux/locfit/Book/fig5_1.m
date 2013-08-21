% Local Regression and Likelihood, Figure 5.1.
% Author: Catherine Loader
%
% Density Estimation, for the Old Faithful Geyser Data.

load geyser;
fit = locfit(geyser,'alpha',[0.1 0.8],'ll',1,'ur',6);
figure('Name','fig5_1: Density estimation for the Old Faithful Geyser Data.' );
lfplot(fit);
xlabel('Old Faithful Eruption Duration');
ylabel('Density');
