% Local Regression and Likelihood, Figure 6.2.
% Author: Catherine Loader
%
% Local smooth of CO2 dataset, main trend.

load co2;
fit = locfit(year+month/12,co2,'alpha',0.5,'deg',1);
figure('Name','fig6_2: CO2 Dataset: Local smoothing, main trend' );
lfplot(fit);
xlabel('Date');
ylabel('CO2')
