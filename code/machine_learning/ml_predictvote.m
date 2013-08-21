function pred = ml_predictvote(trials, model)
% Meta-Prediction function for Voting.
% Prediction = ml_predictvote(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainvote
%
% Out:
%   Prediction  : discrete probability distribution, formatted as
%                 {'disc' [NxC] [Cx1]}, with element #2 being the per-class probability and 
%                 element #3 the original target values per class
%                 thus, the expected target values are Prediction{2}*Prediction{3}
%
% See also:
%   ml_trainvote
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-06-25

if isfield(model,'voted')
    % get the prediction function
    predictor = model.opts{2};

    % get the trial count
    if isnumeric(trials)
        if model.feature_matrices
            trialcount = size(trials,3);
        else
            trialcount = size(trials,1);
        end
    else
        % custom data
        trialcount = utl_default_partitioner(trials);
    end
        
    % construct discrete probability distribution
    pred = {'disc' , zeros(trialcount,length(model.classes)), model.classes};
    
    if strcmp(model.scheme,'1v1')
        % 1v1 voting, adding up the probabilities from each vote
        for i=1:length(model.classes)
            for j=i+1:length(model.classes)
                outcome = predictor(trials,model.voted{i,j});
                if iscell(outcome) && strcmp(outcome{1},'disc')
                    % the classifier produces a discrete probability distribution (preferred)
                    pred{2}(:,[i j]) = pred{2}(:,[i j]) + outcome{2};
                elseif isnumeric(outcome)
                    pred{2}(:,[i j]) = pred{2}(:,[i j]) + [outcome == model.classes(i), outcome == model.classes(j)];
                else
                    error('Unsupported classifier output format');
                end
            end
        end
    else
        % 1vR voting
        classvec = 1:length(model.classes);
        for i=classvec
            % predict i vs rest probabilities
            outcome = predictor(trials,model.voted{i});            
            if iscell(outcome) && strcmp(outcome{1},'disc')
                % the classifier produces a discrete probability distribution (preferred)
                binary_probs = outcome{2};
            elseif isnumeric(outcome)
                % the classifier produces class label estimates
                binary_probs = [outcome == model.classes(i), outcome ~= model.classes(i)];
            else
                error('Unsupported classifier output format');
            end
            % add probabilities for class i
            pred{2}(:,i) = pred{2}(:,i) + binary_probs(:,1);
        end
    end
    
    % renormalize probabilities
    pred{2} = pred{2} ./ repmat(sum(pred{2},2),1,size(pred{2},2));
else
    error('The given model was apparently not constructed via voting.');
end