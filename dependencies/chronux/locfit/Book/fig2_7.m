% Local Regression and Likelihood, Figure 2.7.
%
% Example gcvplot. First argument to gcvplot is a
% vector (or matrix) of smoothing parameters.
% Remaining arguments are passed directly to locfit().
%
% Author: Catherine Loader.
%
% NEED: cpplot.

load ethanol;

alp = (0.2:0.05:0.8)';
figure('Name','fig2_7: Example gcvplot');
gcvplot(alp,E,NOx);
