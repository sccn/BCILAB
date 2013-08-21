function pred = ml_predictsvmperf(trials, model)
% Prediction function for the Support Vector Machine (SVMperf).
% Prediction = ml_predictsvmperf(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainsvmlight
%
% Out:
%   Prediction  : discrete probability distribution, formatted as
%                 {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and 
%                 element #3 the original target values per class
%                 thus, the expected target values are Prediction{2}*Prediction{3}
%
% See also:
%   ml_trainsvmperf
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-04

if isfield(model,'voted')
    pred = ml_predictvote(trials,model);
else
    trials = hlp_applyscaling(trials,model.sc_info);
    raw_pred = min(+1,max(-1,svmperfclassify(trials, ones(size(trials,1),1), model, '-v 0')));
    pred = {'disc', [(1-raw_pred)/2 1-(1-raw_pred)/2], model.classes};
end