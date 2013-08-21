function pred = ml_predictsvm(trials, model)
% Prediction function for the Support Vector Machine.
% Prediction = ml_predictsvmlinear(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainsvmlinear
%
% Out:
%   Prediction  : discrete probability distribution, formatted as
%                 {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and
%                 element #3 the original target values per class
%                 thus, the expected target values are Prediction{2}*Prediction{3}
%
% See also:
%   ml_trainsvm
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-04

if isfield(model,'voted')
    pred = ml_predictvote(trials,model);
else
    % scale the data
    trials = hlp_applyscaling(trials,model.sc_info);
    % kernelize the data
    trials = utl_kernelize(trials,model.basis,model.kernel,model.gammap,model.degree);
    % get raw predicted classes
    if strcmp(model.variant,'native')
        % MATLAB implementation
        class = trials*model.w + model.b > 0;
    else
        % LIBLINEAR implementation
        class = llwpredict(zeros(size(trials,1),1),sparse(double(trials)),model);
    end
    % translate the predicted classes into a discrete probability distribution (for consistency with other classifiers)
    probs = zeros(length(model.classes),length(class));
    probs(length(model.classes)*(0:length(class)-1)'+ class + 1) = 1;
    pred = {'disc', probs', model.classes};
end
