function pred = ml_predictsvmlight(trials, model)
% Prediction function for the Support Vector Machine (SVMlight).
% Prediction = ml_predictsvmlight(Trials, Model)
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
%   ml_trainsvmlight
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-05

if isfield(model,'voted')
    pred = ml_predictvote(trials,model);
else
    trials = hlp_applyscaling(trials,model.sc_info);
    [err,raw_pred] = svmclassify(double(trials), ones(size(trials,1),1), model); %#ok<ASGLU>
    raw_pred = min(+1,max(-1,raw_pred));
    pred = {'disc', [(1-raw_pred)/2 1-(1-raw_pred)/2], model.classes};
end