function pred = ml_predicthkl(trials, model)
% Prediction function for Hierarchical Kernel Learning.
% Prediction = ml_predicthkl(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainhkl
%
% Out:
%   Prediction  : in the case of classification ('bernoulli' likelihood was used):
%                   discrete probability distribution, formatted as
%                   {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and 
%                   element #3 the original target values per class
%                   thus, the expected target values are Prediction{2}*Prediction{3}
%                 in the case of regression ('gaussian' or other likelihood was used):
%                   [Nx1] vector of predicted target values%
% See also:
%   ml_trainhkl
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-04

if isfield(model,'voted')
    pred = ml_predictvote(trials,model);
else
    % predict
    stats = hkl_test(model.model,model.outputs,trials,'Ytest',ones(size(trials,1),1));
    pred = stats.predtest(:,model.bestlambda);
        
    % in case of classification: map to probabilities and create a discrete distribution
    if strcmp(model.loss,'logistic')
        pred = logsig(pred-0.5);
        pred = {'disc' [1-pred pred] model.classes};
    end
end