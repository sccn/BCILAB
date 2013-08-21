% Local Regression and Likelihood, Figure 7.2.
%
% Conditional Hazard Rate Estimation.
% When hazard rate estimation is performed with multiple x variables,
% the result is an estimate of the conditional hazard rate for the
% first variable, given the levels of the remaining variables.
%
% The surface plot looks more dramatic than the contour plot in the book!
%
% Author: Catherine Loader

load livmet;
fit = locfit([t dm],'cens',1-z,'scale',0,'deg',1,'family','hazard','alpha',0.5,'xlim',[0 1;10000 20]);
figure('Name','fig7_2: Conditional hazard rate estimation' );
lfplot(fit);
xlabel('Survival Time (Months)');
ylabel('Diameter (c.m.)');
zlabel('Hazard Rate');
