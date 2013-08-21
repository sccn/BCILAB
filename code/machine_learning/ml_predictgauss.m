function pred = ml_predictgauss(trials, model)
% Prediction function for the Gaussian Bayes classifier.
% Prediction = ml_predictgauss(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_traingauss
%
% Out:
%   Prediction  : discrete probability distribution, formatted as
%                 {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and 
%                 element #3 the original target values per class; probabilities are unnormalized.
%
% See also:
%   ml_traingauss
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-05

% scale data & and select features
trials = hlp_applyscaling(trials,model.sc_info);

if model.normprobs >= 0
    trials = trials(:,model.feature_mask);

    % compute raw per-class probabilities...
    probs = zeros(size(trials,1),length(model.classes));
    for c=1:length(model.classes)
        probs(:,c) = mvnpdf(trials,model.mu{c},model.sigma{c}); end
    
    % optionally renormalize them
    if model.normprobs
        probs = probs ./ repmat(sum(probs,2),1,size(probs,2)); end
    
    pred = {'disc', probs, model.classes};
else
    % (undocumented feature: use the feature-selection model directly...)
    pred = ml_predictlogreg(trials,model.selector); 
end
