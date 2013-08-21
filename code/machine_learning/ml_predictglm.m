function pred = ml_predictglm(trials,model)
% Simple prediction function for the Bayesian GLM.
% Prediction = ml_predictglm(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_traindal
%
% Out:
%   Prediction  : discrete probability distribution, formatted as
%                 {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and 
%                 element #3 the original target values per class
%                 thus, the expected target values are Prediction{2}*Prediction{3}
%
% See also:
%   ml_trainglm
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-25

% vectorize data if necessary
if model.vectorize
    trials = double(reshape(trials,[],size(trials,3))'); end

if isfield(model,'voted')
    % dispatch to the voter
    pred = ml_predictvote(trials,model);
else
    % scale data
    trials = hlp_applyscaling(trials,model.sc_info);
    
    % add bias term
    if isfield(model,'b') && ~isempty(model.b)
        model.w(end+1) = model.b;    
        trials = [trials ones(size(trials,1),1)];
    end
    
    if strcmp(model.ptype,'classification')
        % map to probabilities
        probs = 1 ./ (1 + exp(-trials*full(model.w)));
        pred = {'disc', [1-probs probs], model.classes};
    else
        % do linear regression
        pred = trials*full(model.w);
    end
end
