% Model the success probability of successive trials of a monkey
% performing a task.
%
% The 'y' variable is a vector of 0/1, with 1 denoting success on
% the trial; 0 failure. The fitting procedure uses logistic
% regression within sliding windows (specified by the 'family','binomial'
% arguments).
%
% Choosing the bandwidth here is critical. The data shows `on/off' behavior,
% exhibiting periods of mainly successes, and mainly failures, respectively.
% Large values of alpha will smooth out this behavior, while small values
% will be too sensitive to random variability. Values of 0.15 to 0.2 seem
% reasonable for this example.
%
% AIC is based on asymptotic approximations, and seems unreliable here --
% formal model selection needs more investigation.
%
% Data is from Keith Purpura.

load 050527_correct.mat;
y = byTrial(1).correct';
n = length(y);
fit = locfit((1:n)',y,'family','binomial','alpha',0.15);
lfplot(fit);
title('Local Logistic Regression - Estimating Success Probability');
lfband(fit);