% Local Regression and Likelihood, Figure 9.2.
%
% Hardle's Motorcycle accelaration dataset.
% Local Variance Estimation.
%
% Author: Catherine Loader
%
% NEED: lfknots function to extract fit points.

load mcyc;
fit = locfit(time,'y',accel,'nn',0.1);
y = -2*predict(fit,'fitp','what','lik');
w = predict(fit,'fitp','what','rdf');
x = fit.fit_points.evaluation_points';
[x y w]
fitv = locfit(x,y,'weights',w,'family','gamma','nn',0.4);
figure('Name','fig9_2: Motorcycle acceleration');
lfplot(fitv);
xlabel('Time');
ylabel('Variance');
