function pred = ml_predictdal(trials,model)
% Prediction function for Dual-Augmented Lagrangian.
% Prediction = ml_predictdal(Trials, Model)
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
%   ml_traindal
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-25

if isfield(model,'voted')
    % dispatch to the voter
    pred = ml_predictvote(trials,model);
else
    % vectorize data if necessary
    if model.vectorize
        trials = double(reshape(trials,[],size(trials,3))'); end

    % scale data
    trials = hlp_applyscaling(trials,model.sc_info);

    % transpose if necessary
    if model.transpose
        trials = double(reshape(trials',model.shape(2),model.shape(1),[]));
        ntrials = zeros(model.shape(1),model.shape(2),size(trials,3));
        for t=1:size(trials,3)
            ntrials(:,:,t) = trials(:,:,t)'; end
        trials = double(reshape(ntrials,[],size(ntrials,3))');
    end
    
    % add bias term to data
    trials = [trials ones(size(trials,1),1)];

    % add bias to model parameters, too
    model.w = model.w(:);
    model.w(end+1) = model.b;
    
    if strcmp(model.loss,'logistic')
        % map to probabilities
        probs = 1 ./ (1 + exp(-trials*full(model.w)));
        pred = {'disc', [1-probs probs], model.classes};
    else
        % do linear regression
        pred = trials*full(model.w);
    end
end