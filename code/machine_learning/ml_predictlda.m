function pred = ml_predictlda(trials, model)
% Prediction function for Linear Discriminant Analysis
% Prediction = ml_predictlda(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainlda
%
% Out:
%   Prediction  : discrete probability distribution, formatted as
%                 {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and 
%                 element #3 the original target values per class
%                 thus, the expected target values are Prediction{2}*Prediction{3}
%
% Examples:
%   targets might look like this: [-1 -1 1 -1 1 -1 -1 1 -1 -1 1 -1 -1 1 ...]' 
%
%   model = ml_trainlda(data,targets)
%   p = ml_predictlda(data, model); expectation = p{2}*p{3};
%   now expectation might look like this: [-0.6 -0.9 +0.4 -0.7 +0.8 -0.1 +0.5 +1.0 -0.9 +1.0 -1.0 -1.0 +1.0 ...]'
%
% See also:
%   ml_trainlda
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-03


if isfield(model,'voted')
    % dispatch to the voter for multi-class classification
    pred = ml_predictvote(trials,model);
else
    % otherwise calculate the labels
    trials = trials(:,model.featuremask);
    raw_labels = min(+1,max(-1,trials*model.w' - model.b));
    pred = {'disc', [(1-raw_labels)/2 1-(1-raw_labels)/2], model.classes};
end