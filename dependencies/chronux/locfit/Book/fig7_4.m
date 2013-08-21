% Local Regression and Likelihood, Figure 7.4.
%
% Censored Local Likelihood
%
% Author: Catherine Loader
%
% NEED: this is given problems for unknown reasons. something
% to do with censored observations and the ibeta routine.
% if i replace this by
% res[ZLIK] = th*y-y*log(p)
% res[ZDLL] = y*p
% res[ZDDLL] = y*p*(1-p)
% works fine. Memory leak somewhere??
% also works if i ignore censoring.


load border;
fit = locfit(day,runs,'cens',no,'family','geom','alpha',0.7);
figure('Name','fig7_4: censored local likelihood;' );
lfplot(fit);
xlabel('Date');
ylabel('Runs');
