function pred = ml_predictproximal(trials,model)
% Prediction function for the proximal framework.
% Prediction = ml_predictproximal(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainproximal
%
% Out:
%   Prediction  : discrete probability distribution, formatted as
%                 {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and 
%                 element #3 the original target values per class
%                 thus, the expected target values are Prediction{2}*Prediction{3}
%
% See also:
%   ml_trainproximal
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-02-05

if isfield(model,'voted')
    % dispatch to the voter
    pred = ml_predictvote(trials,model);
else
    % vectorize data if necessary
    if model.vectorize_trials
        trials = double(reshape(trials,prod(model.featureshape),[])'); end

    % scale data
    trials = hlp_applyscaling(trials,model.sc_info);

    % add bias term to data
    if model.includebias
        trials = [trials ones(size(trials,1),1)]; end

    model.w = model.w(:);
    
    if strcmp(model.loss,'logistic')
        % map to probabilities
        probs = 1 ./ (1 + exp(-trials*full(model.w)));
        if ~model.continuous_targets
            pred = {'disc', [1-probs probs], model.classes};
        else
            pred = probs;
        end
    else
        % do linear regression
        pred = trials*full(model.w);
    end
end
