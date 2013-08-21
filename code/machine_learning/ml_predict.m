function pred = ml_predict(trials, model)
% Make predictions for some data, using some (previously learned) model.
% Prediction = ml_predict(Trials, Model)
%
% This function makes predictions for some given data and a supplied model, as learned using
% ml_train. It dispatches to the appropriate ml_predict* function, depending on what was used to
% train the model.
%
% For each of the N supplied trials, one prediction is produced in the output, which may be in a
% number of different formats, depending on the parameters of the model (which in turn depends on
% the learning function used, its parameters, and the form of the target variables in the training
% data).
%
% The most common format of the prediction is discrete probability distributions per trial, which is
% the output produced by most classifiers. In this case, structed as a 3-element cell array {'disc',
% Probabilities, Classes}, the distributions are returned as a [#Trials x #Classes] matrix in the
% Probabilities entry, which contains the probability for each of a set of classes, per trial. The
% set of classes is associated with set of target values (e.g. -1,+1, or 3,7,5) by the Classes cell
% entry, which is size [#Classes x 1], and is sorted in the same order as in the probability
% distributions (in ascending order of target values), so that Prediction{2}*Prediction{3} would
% give a probability-weighted sum of the possible target values). This format can be mapped into the
% "usual" format of most likely target value per trial via the expression
% Prediction{3}(argmax(Prediction{2}')).
%
% For regression, the most common format is the [NxD] (usually D=1) format of one point estimate per
% trial, since only very few methods can give a posterior probability distribution right now. Other
% formats are primarily provided for more advanced predictors which can give regression outputs as a
% parametric (or even non-parametric) probability distribution per trial, as {distrib,NxP}, or in an
% entirely custom format, as {'struct',{N}}.
%
% In:
%   Trials   : data, [NxF] array (N = number of training instances,F = number of feature dimensions) 
%              or in special cases a [UxVxN] array (N... number of training instances, U,V ...
%              number of rows/columns of the feature matrices)
%              If not otherwise possible, it is also allowed to pass trials in any custom format to
%              ml_predict, given that the prediction function in question can handle these data.
%
%   Model    : the model to be used for prediction (output of ml_train);
%              when empty, the actual features are returned (passed through).
%
% Out:
%   Prediction : predictions of the model, one for each trial
%                format depends on the chosen Model:
%                  * [Nx1] : classification or one-dimensional regression outputs (point estimates)
%                  * [NxD] : d-dimensional regression outputs (point estimates)
%                  * {'disc', [NxC], [Cx1]}: discrete probabilities for C classes in probabilistic
%                    classification, the last array determines the assigment from class indices (in
%                    NxC) to the output labels, such that [NxC]*[Cx1] is the expected output
%                  * {distrib, NxP}: posterior probability distributions in probabilistic
%                    regression, with distrib being one of:
%                    'bino','chi2','exp','ev','f','gam','gev','gp','geo','hyge','logn','nbin','ncf',
%                    'nct','ncx2','norm','poiss','rayl','t','unif','unid','wbl' with the
%                    appropriate interpretation under pdf(), and P the number of parameters to pdf()
%                  * {distrib, {N}}: posterior probability distributions, with distrib being either
%                    'mvn', for mvnpdf(), or any function handle that can accept parameter sets of
%                    the form distrib(X,Prediction{2}{i}{:});
%                  * {'struct', {N}}: structured predictions
%
% See also:
%   ml_train
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-03

if isempty(model)
    pred = trials;
else   
    pred = feval(['ml_predict' model.args.arg_selection],trials,model.model);
    
    % post-process predictions
    if iscell(pred) && length(pred) == 3 && strcmp(pred{1},'disc')
        % sort discrete distributions
        [classes,inds] = sort(pred{3});
        pred = {pred{1} pred{2}(:,inds) classes};
    end
end