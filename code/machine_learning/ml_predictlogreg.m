function pred = ml_predictlogreg(trials,model)
% Prediction function for Logistic Regression.
% Prediction = ml_predictlogreg(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainlogreg
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
%   model = ml_trainlogreg(data,targets)
%   p = ml_predictlogreg(data, model); expectation = p{2}*p{3};
%   now expectation might look like this: [-0.6 -0.9 +0.4 -0.7 +0.8 -0.1 +0.5 +1.0 -0.9 +1.0 -1.0 -1.0 +1.0 ...]'
%
% See also:
%   ml_trainlogreg
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-04

if isfield(model,'voted')
    % dispatch to the voter
    pred = ml_predictvote(trials,model);
else    
    % scale data
    trials = hlp_applyscaling(trials,model.sc_info);
    
    % add bias variable
    if ~isempty(strfind(model.variant,'vb'))
        trials = [trials ones(size(trials,1),1)]; end
    
    switch model.variant
        case {'vb','vb-iter','vb-ard'}
            probs = bayes_logit_post(trials, model.w, model.V, model.invV);
            if ~model.continuous_targets
                pred = {'disc', [1-probs probs], model.classes};
            else
                pred = probs;
            end
        case 'lars'
            probs = real(glmnetPredict(model.CVerr.glmnet_object,'response',double(trials),model.CVerr.lambda_min));
            if ~model.continuous_targets
                pred = {'disc', probs, model.classes};
            else
                pred = probs;
            end
        case 'l1'
            if isfield(model,{'W','b'})
                probs = 1 ./ (1 + exp(-trials*model.W + model.b));
                pred = {'disc', [1-probs probs], model.classes};
            else
                if model.bias
                    probs = 1 ./ (1 + exp(-trials*model.w(1:end-1)' - model.w(end)));
                else
                    probs = 1 ./ (1 + exp(-trials*model.w'));
                end
                % doesn't give probabilistic outputs: [x,y,probs] = llpredict(zeros(size(trials,1),1),sparse(double(trials)),model,'-b 1'); %#ok<ASGLU>
                pred = {'disc', [probs 1-probs], model.classes(model.Label+1)};
            end
        case 'l2'
            if isfield(model,{'W','b'})
                probs = 1 ./ (1 + exp(-trials*model.W + model.b));
                pred = {'disc', [1-probs probs], model.classes};
            else
                [x,y,probs] = llwpredict(zeros(size(trials,1),1),sparse(double(trials)),model,'-b 1'); %#ok<ASGLU>
                pred = {'disc', probs, model.classes(model.Label+1)};
            end
    end
end