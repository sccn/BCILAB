function pred = ml_predictrvm(trials,model)
% Prediction function for the Relevance Vector Machine.
% Prediction = ml_predictrvm(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainrvm
%
% Out:
%   Prediction  : in the case of classification ('bernoulli' likelihood was used):
%                   discrete probability distribution, formatted as
%                   {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and 
%                   element #3 the original target values per class
%                   thus, the expected target values are Prediction{2}*Prediction{3}
%                 in the case of regression ('gaussian' or other likelihood was used):
%                   [Nx1] vector of predicted target values
%
% See also:
%   ml_trainrvm
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-06

if isfield(model,'voted')
    pred = ml_predictvote(trials,model);
else    
    % scale the data
    trials = hlp_applyscaling(trials,model.sc_info);
    
    % kernelize it
    trials = utl_kernelize(trials,model.basis,model.kernel,model.gamma,model.degree);
    
    % add a bias if necessary
    if model.bias        
        trials = [ones(size(trials,1),1) trials]; end
    
    if strcmp(model.kernel,'linear')
        trials = trials(:,model.feature_sel); end
    
    % map to predictions
    pred = trials*model.param.Value;
    
    % send through the link function if necessary
    if model.ptype == 'c'
        pred = SB2_Sigmoid(pred);
        pred = {'disc' [1-pred pred] model.classes};
    end
end