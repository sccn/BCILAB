function model = ml_traincov(varargin)
% Learn a linear predictive model using covariance-based classification.
% Model = ml_traincov(Trials, Targets, Options...)
%
% This method assumes that the given trials values represent covariance matrices and offers various
% methods, including information geometric approaches, to classify them. This implementation uses
% the covariance toolbox by Alexandre Barachant.
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%
% Out:
%   Model   : the predictive model; can be used with ml_predictcov
%
% Examples:
%   % learn a model using the defaults
%   model = ml_traincov(trials,targets);
%
% See also:
%   ml_predictcov
%
%                           Christian Kothe, Syntrogi
%                           2015-04-21
        
arg_define([0 3],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'classifier','ClassificationMethod'},'mdm',{'mdm','fgmdm','tslda'}, 'Classification method to use. The mdm method (minimum distance to mean) is based on minimum distance to class mean using a given distance metric and mean estimation metric. The fgmdm (filtered geodesic mdm) method additionally performs geodesic filtering. The tslda method (tangent-space linear discriminant analysis) is an LDA implementation that respects the Riemannian geometry of the data.'), ...
    arg({'mean_est','MeanEstimator'},'riemann',{'arithmetic','riemann','riemanndiag','riemanntrim','median','riemannmed','logeuclid','opttransp','ld','geodesic','harmonic','geometric'},'Method to average covariance matrices. Various methods are supported, including Riemannian mean/median/trimmed mean, log-euclidean mean, optimal transportation mean, log determinant mean, and geodesic iterative mean.'), ...
    arg({'distance_metric','DistanceMetric'},'riemann',{'euclid','riemann','kullback','logeuclid','opttransp','ld'},'Distance metric to use. This is only used for the mdm and fgmdm methods. Different distance metrics are supported, including the Euclidean metric, the Riemannian distance, Kullback-Leibler divergence, log-euclidean distance, optimal transportation distance, and log determinant distance'));

% find the class labels
classes = unique(targets);
if length(classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
else
    % reformat feature shape to DxDxN
    if ndims(trials) == 2 %#ok<ISMAT,NODEF>
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
    
    model.classes = classes;
    switch classifier
        case 'mdm'
            % minimum distance to mean: estimate class means
            for c=length(classes):-1:1
                model.C{c} = mean_covariances(trials(:,:,targets==classes(c)),mean_est); end
        case 'fgmdm'
            % filtered geodesic minimum distance to mean                        
            % geodesic filtering
            [model.W,model.Cg] = fgda(trials,targets,mean_est,{},'shcov',{});
            trials = geodesic_filter(trials,model.Cg,model.W(:,1:length(model.classes)-1));
            % estimate class means
            for c=length(classes):-1:1
                model.C{c} = mean_covariances(trials(:,:,targets==classes(c)),mean_est); end
        case 'tslda'
            error('This variant is not yet implemented!');
        otherwise
            error('Unsupported classification method: %s',hlp_tostring(classifier));
    end
end

model.classifier = classifier;
model.mean_est = mean_est;
model.distance_metric = distance_metric;
