function pred = ml_predictqda(trials, model)
% Prediction function for Quadratic Discriminant Analysis.
% Prediction = ml_predictqda(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainqda
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
%   model = ml_trainqda(data,targets)
%   p = ml_predictqda(data, model); expectation = p{2}*p{3};
%   now expectation might look like this: [-0.6 -0.9 +0.4 -0.7 +0.8 -0.1 +0.5 +1.0 -0.9 +1.0 -1.0 -1.0 +1.0 ...]'
%
% See also:
%   ml_trainqda
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-03

if isfield(model,'voted')
    % dispatch to the voter
    pred = ml_predictvote(trials,model);
else
    % pre-prune features
    trials = trials(:,model.featuremask);
    % pre-scale the data
    trials = hlp_applyscaling(trials,model.sc_info);    
    % calculate the labels
    raw_labels = min(+1,max(-1,-(sum(((trials*model.q).*trials)') + model.l*trials' - model.c)))'; %#ok<UDIM>
    pred = {'disc', [(1-raw_labels)/2 1-(1-raw_labels)/2], model.classes};
end