function model = ml_trainhkl(varargin)
% Learn a sparse non-linear predictive model using Hierarchical Kernel Learning.
% Model = ml_trainhkl(Trials, Targets, Lambdas, Options...)
%
% Hiararchical Kernel Learning [1,2] is a non-linear, kernel-based (as in Support Vector Machines)
% learning algorithm which can be used both for regression and classification. It can give
% probabilistic outputs (similar to logistic regression). The defining property of HKL is that it
% can select relevant features out of a large number of potentially many irrelevant ones (similar to
% LASSO), whereas these features are selected according to non-linear criteria, which makes HKL
% uniquely powerful in selecting relevant portions of the data.
%
% The method is relatively slow due to the need for regularization, but on the other hand has the
% highest chance of learning relevant structure in very large numbers of dimensions (among the
% provided machine learning methods). For these reasons, it is not the best method to start with, as
% others allow to evaluate data and to iterate through parameters at a much faster rate. It is,
% however, one of the best methods to try in order to get the best possible results, and also to
% check whether data on which other methods have failed contains any information of interest at all.
% 
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%
%   Lambdas      : vector of regularization parameters to optimize for
%                  default: 10.^[1:-.5:-8]; note: if this is specified as a vector instead of a search() expression,
%                  ml_trainhkl will internally optimize the parameter using a more efficient procedure; search() can nevertheless
%                  be used instead, to get finer control over cross-validation (e.g. safety margins, etc.)
%
%   Options  : optional name-value pairs to control the training details:
%              'loss'  : either 'logistic' for classification or 'squared' for regression (default: 'logistic')
%              'kernel': type of kernel to use; can be one of the following;
%                        'hermite': hermite polynomial kernel (default)
%                        'spline' spline kernel 
%                        'gauss-hermite-full':  hermite expansion of the Gaussian kernel (finite order)
%                        'anova': ANOVA kernel
%                        'polynomial': polynomial kernel (not recommended)
%              'params': vector of kernel parameters; if empty, chosen according to the author's recommendations (default: [])
%                        see hkl.m for a detailed explanation of this parameter
%              'mingap': result precision; expressed as duality gap (default: 1e-3)
%              'gapeta': smoothing of the reduced problem (default: 1e-3)
%              'maxactive': maximum number of selected kernel (default: 100)
%              'display': display diagnostic outputs (default: 0)
%              'data_normalization': either 'center' (for centering) or 'scale' (for centering and standardization) (default: 'scale')
%              'kernel_normalization': 'center' (for centering) 'scale' (for centering and standardization), or 'scale-root' (default: 'scale-root')
%              'alpha': dual parameter, see [2]
%              'eta': kernel weights, see [2]
%              'b': constant term, see [2]
%              'memory_cache': size of the memory cache; should be large enough to fill most of the free RAM (default: 1e9)
%              'k': foldness of an internal cross-validation, to compute the optimal regularization parameter (default: 5)
%
% Out:
%   Model   : a HKL model...
%
% See also:
%   ml_predicthkl 
%
% Examples:
%   % learn a standard HKL classifier
%   model = ml_trainhkl(trials,targets);
%
%   % as before, but use a specific set of regularization values
%   model = ml_trainhkl(trials,targets.10.^[1:-1:-8]);
%
%   % like before, but use the squared loss for regression
%   model = ml_trainhkl(trials,targets.10.^[1:-1:-8],'loss','squared');
%
%   % like before, but this time use a different kernel type
%   model = ml_trainhkl(trials,targets.10.^[1:-1:-8],'kernel','spline');
%
%
% References:
%   [1] F. Bach. "Exploring Large Feature Spaces with Hierarchical Multiple Kernel Learning." 
%       Advances in Neural Information Processing Systems (NIPS), 2008.
%   [2] F. Bach. "High-Dimensional Non-Linear Variable Selection through Hierarchical Kernel Learning." 
%       Technical report, HAL 00413473, 2009.
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-04

opts=arg_define([0 3],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'lambdas','Lambdas'}, 10.^(1:-.5:-8), [0 10^2 10^-10 Inf], 'Vector of regularization parameters to consider. If this is specified as a vector instead of a search() expression, HKL will internally optimize the parameter using a more efficient procedure than the framework''s search() facility.','cat','Core Parameters'), ...
    arg({'loss','LossFunction'}, 'logistic', {'logistic','squared'}, 'Loss function to be used. The logistic loss is suited for classification problems, whereas the squared loss is suited for regression problems.','cat','Core Parameters'), ...
    arg({'kernel','KernelFunc'}, 'hermite', {'hermite','spline','gauss-hermite-full','anova','polynomial'}, 'Class of kernels to use. Hermite polynomials, splines, Hermite expansion of the Gaussian kernel, ANOVA kernel, or the (not recommended) polynomial kernel.','cat','Core Parameters'), ...
    arg({'params','KernelParams'}, [], [], 'Vector of kernel parameters. If empty, chosen according to the author''s recommendations (Francis Bach).','shape','row','cat','Core Parameters','guru',true), ...
    arg({'mingap','Precision'}, 1e-3, [], 'Desired precision. Expressed as the duality gap.','cat','Core Parameters','guru',true), ...
    arg({'gapeta','GapETA'},1e-3,[],'Smoothing degree of the reduced problem.','cat','Core Parameters','guru',true), ...
    ...
    arg({'maxactive','MaxKernelsActive'},400,[0 Inf],'Maximum active kernels..','cat','Efficiency'), ...
    arg({'memory_cache','MemoryCache'},0.25,[0.01 Inf],'Size of the memory cache. If < 1, assumed to be a fraction of the current free memory, otherwise in bytes. Should be large enough to fill most of the free RAM for efficiency.','cat','Efficiency'), ...
    arg({'k','Foldness'},5,[0 Inf],'Cross-validation folds for parameter search.','cat','Efficiency'), ...
    ...
    arg({'showit','Verbose','display'},true,[],'Display diagnostic output.','cat','Miscellaneous'), ...
    arg({'data_normaliation','DataNormalization'},'center',{'center','scale'},'Data normalization. Very important if the data is not naturally normalized. Note: scale both centers and standardizes.','cat','Miscellaneous','guru',true), ...
    arg({'kernel_normaliation','KernelNormalization'},'scale-root',{'center','scale','scale-root'},'Normalization in the kernel space.','cat','Miscellaneous','guru',true), ...
    arg({'alphaparam','Alpha','alpha'},[],[],'Dual parameter. See the HAL 00413473 TechReport.','cat','Miscellaneous','guru',true), ...
    arg({'eta','Eta'},[],[],'Kernel Weights. See the HAL 00413473 TechReport.','cat','Miscellaneous','guru',true), ...
    arg({'b','B'},[],[],'Constant Term. See the HAL 00413473 TechReport.','cat','Miscellaneous'));

arg_toworkspace(opts);

if isempty(lambdas)
    lambdas = 10.^(1:-.5:-8); end

if iscell(targets) %#ok<*NODEF>
    targets = targets{1};
    disp('note: hkl does not support weighted learning.');
end

% identify and remap the classes, if necessary
classes = unique(targets);
if strcmp(loss,'logistic')
    if length(classes) > 2
        % in this case we use the voter...
        model = ml_trainvote(trials,targets,'1v1',@ml_trainhkl,@ml_predicthkl,varargin{:});
        return;
    elseif length(classes) == 1
        error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
    else
        % remap targets
        targets(targets == classes(1)) = 0;
        targets(targets == classes(2)) = 1;
    end
end

% select default kernel parameter, if necessary
params = hlp_rewrite(kernel,'polynomial',[3 .1 4], 'hermite',[.5 3 .1 4], 'gauss-hermite',[1 .05 3 .1 .5], 'gauss-hermite-full',[1 .05 3 .1 .5 30], 'anova',[.05 .1 4 30], 'spline',[.1 4 40]);


% turn into fraction of free memory
if opts.memory_cache < 1
    opts.memory_cache = hlp_memfree * opts.memory_cache; end

% learn
warning off MATLAB:illConditionedMatrix
args = hlp_struct2varargin(opts,'rewrite',{'alphaparam','alpha'},'restrict', {'mingap','gapeta','maxactive','data_normalization','kernel_normalization','alpha','eta','b','memory_cache','conjgrad'});
if length(lambdas)>1
    % find the best lambda by randomized CV
    if ~showit
        evalc('[model.outputs,model.model,dum1,dum2,model.bestlambda] = hkl_kfold(k,''same'',trials,targets,lambdas,loss,kernel,params,args{:});');
    else
        [model.outputs,model.model,dum1,dum2,model.bestlambda] = hkl_kfold(k,'same',trials,targets,lambdas,loss,kernel,params,args{:}); %#ok<ASGLU>
    end
else
    % compute the model directly for the given lambda
    [model.outputs,model.model] = hlp_diskcache('predictivemodels',@hkl,trials,targets,lambdas,loss,kernel,params,args{:});
    model.bestlambda = 1;
end

model.classes = classes;
model.loss = loss;
