% Local Regression and Likelihood, Figure 7.1.
%
% Estimating the hazard rate for censored survival (or
% failure time) data. Hazard rate estimation is specified
% by 'family','hazard'.
%
% The censoring indicator variable is passed as the 'cens'
% argument. This should be a vector of 0's and 1's, with 1
% indicating a censored observation, 0 uncensored.
%
% The 'xlim' argument specifies bounds on the domain of the
% x variable (survival times). Usually, survival times are
% non-negative (lower bound 0), but with theoretically no
% upper bound. The specification [0;10000] gives the lower
% bound of 0, and the upper bound is effectively infinite.
%
% Author: Catherine Loader
%
% NEEDS: m argument to lfplot.
% Also, lfplot symbols should distinguish between censored and
% uncensored data points.

load heart;
fit = locfit(surv,'cens',cens,'family','hazard','alpha',0.4,'xlim',[0;10000]);
figure('Name','fig7_1: Hazard rate for censored data' );
lfplot(fit);
xlabel('Survival Time');
ylabel('Hazard Rate');
