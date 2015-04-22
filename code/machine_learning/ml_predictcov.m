function pred = ml_predictcov(trials, model)
% Prediction function for covariance-based classification
% Prediction = ml_predictcov(Trials, Model)
%
% In:
%   Trials  : the data a matrix, as in ml_predict
%
%   Model   : predictive model as produced by ml_traincov
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
%   model = ml_traincov(data,targets)
%   p = ml_predictcov(data, model); expectation = p{2}*p{3};
%   now expectation might look like this: [-0.6 -0.9 +0.4 -0.7 +0.8 -0.1 +0.5 +1.0 -0.9 +1.0 -1.0 -1.0 +1.0 ...]'
%
% See also:
%   ml_traincov
%
%                           Christian Kothe, Syntrogi
%                           2015-04-21

% reformat feature shape to DxDxN
if ndims(trials) == 2 %#ok<NODEF>
    [N,F] = size(trials);
    D = sqrt(F);
    if abs(D-round(D)) > 0
        error('The number of features in your trials must be a square of an integer.'); end
    trials = reshape(trials',[D,D,N]);
else
    [U,V,N] = size(trials); %#ok<ASGLU>
    if U ~= V
        error('Your feature matrices are not square, i.e., cannot be covariance matrices.'); end
end

% optional geodesic pre-filtering
if strcmp(model.classifier,'fgmdm')
    trials = geodesic_filter(trials,model.Cg,model.W(:,1:length(model.classes)-1)); end

% calculate distances between trials and class centers
Nclass = length(model.classes);
Ntrials = size(trials,3);
d = zeros(Ntrials,Nclass);
for j=1:Ntrials
    for i=1:Nclass
        d(j,i) = distance(trials(:,:,j),model.C{i},model.distance_metric); end
end

% calculate avg distance from a class to all other classes
p = zeros(1,Nclass);
for i=1:Nclass
    for j=Nclass:-1:1
        if j == i
            x(j) = NaN;
        else
            x(j) = distance(model.C{i},model.C{j},model.distance_metric); 
        end
    end
    p(i) = sum(x(~isnan(x)))/(Nclass-1);
end

% normalize distances by that
d = bsxfun(@times,d,1./p);

% calculate class scores between 0 and 1
pot = exp(-d*3); 

% convert to pseudo-probabilities
probs = bsxfun(@times,pot,1./sum(pot,2));

% format predictions
pred = {'disc', probs, model.classes};
