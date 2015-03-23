function pred = ml_predicthdca(trials, model)
% Prediction function for Hierarchical Discriminant Component Analysis
% Prediction = ml_predicthdca(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_trainhdca
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
%   model = ml_trainhdca(data,targets)
%   p = ml_predicthdca(data, model); expectation = p{2}*p{3};
%   now expectation might look like this: [-0.6 -0.9 +0.4 -0.7 +0.8 -0.1 +0.5 +1.0 -0.9 +1.0 -1.0 -1.0 +1.0 ...]'
%
% See also:
%   ml_trainhdca
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-03


if isfield(model,'voted')
    % dispatch to the voter for multi-class classification
    pred = ml_predictvote(trials,model);
else
    % for each modality range...
    for m=1:length(model.modality_ranges)
        range = model.modality_ranges{m};
        
        % for each sub-block...
        for b=size(trials,2):-1:1
            % extract per-block features
            blocktrials = reshape(trials(range,b,:),[],size(trials,3))';
            blockpredictions{b} = blocktrials*model.blockmodels{m}{b}.w' - model.blockmodels{m}{b}.b';
        end
        % extract per-modality features
        layertrials = cat(2,blockpredictions{:});
        rangepredictions{m} = layertrials*model.rangemodels{m}.w' - model.rangemodels{m}.b'; %#ok<AGROW>
    end
    % perform top-level prediction
    toptrials = cat(2,rangepredictions{:});
    pred = ml_predictlda(toptrials,model.topmodel);
end
