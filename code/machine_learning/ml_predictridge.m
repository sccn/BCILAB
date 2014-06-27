function pred = ml_predictridge(trials, model)
% Prediction function for Ridge Regression.
% Prediction = ml_predictridge(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainridge
%
% Out:
%   Prediction  : [Nx1] vector of predicted target values
%
% See also:
%   ml_trainridge
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2013-11-17

trials = hlp_applyscaling(trials,model.sc_info);
trials = [trials ones(size(trials,1),1)];
pred = trials*model.w;
