% Local Regression and Likelihood, Figure 11.1
%
% Variable degree fit, uses module 'vord'. Note there is
% no `lower' degree here; it defaults to 0.
%
% Author: Catherine Loader

load ethanol;

figure('Name','fig11_1a: Variable degree fit');
fit = locfit(E,NOx,'deg',3,'nn',0.3,'module','vord');
lfplot(fit);

figure('Name','fig11_1b: Variable degree fit');
x = fit.fit_points.evaluation_points';
z = predict(fit,'fitp','what','deg');
plot(x,z,'o');
xlabel('Fitting Point');
ylabel('Degree');
