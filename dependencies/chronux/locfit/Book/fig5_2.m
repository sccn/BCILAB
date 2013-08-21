% Local Regression and Likelihood, Figure 5.2.
% Author: Catherine Loader
%
% Density Estimation, for the Old Faithful Geyser Data.
% identity link.

load geyser;
fit = locfit(geyser,'alpha',[0.1 0.8],'link','ident','ll',1,'ur',6);
figure('Name','fig5_2: Density estimation for Old Faithful geyser data' );
lfplot(fit);
xlabel('Old Faithful Eruption Duration');
ylabel('Density');
