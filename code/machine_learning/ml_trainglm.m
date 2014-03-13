function model = ml_trainglm(varargin)
% Learn a Bayesian generalized linear model
% Model = ml_trainglm(Trials, Targets, Options...)
%
% This is an incomplete experimental implementatation wrapping the glm-ie toolbox; more to come
% later.
%
% In:
%   Trials       : training data, as in ml_train
%                  in addition, it may be specified as UxVxN 3d matrix,
%                  with UxV-formatted feature matrices per trial (N trials), or
%                  as {{U1xV1,U2xV2,...}, {U1xV1,U2xV2,...}
%
%   Targets      : target variable, as in ml_train
%
%   Options  : optional name-value parameters to control the training details:
%              'ptype': problem type: 'classification' (default) or 'regression'
%
%              'shape': if trials is a NxF 2d matrix of vectorized matrices,
%                           this is the dimensions of the matrices (default: Fx1)
%                       if trials is specified as an UxVxN 3d matrix, shape defaults to
%                           [U,V] and trials are vectorized into the regular [N, U*V]
%                       if shape is specified with one of the values being NaN,
%                           that value is set so that prod(shape)==F==U*V
%
%              'scaling': pre-scaling of the data (see hlp_findscaling for options) (default: 'std')
%
%              'setupfcn': Setup Function. Receives trials, targets and the options structure, and 
%                          returns the variables [X,y,B,pot,tau,G]. If empty, performs ridge
%                          regression or logistic regression depending on Type parameter.
%
%               for additional parameters see infEngine.m in the glm-ie directory.
%
% Out:
%   Models   : a predictive model
%
% See also:
%   ml_predictglm, dli
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-07-08


args = arg_define([0 2],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    ... % (pre)processing options
    arg({'lambdas','Lambdas'}, 1, [0 Inf], 'Noise variance parameter. This can be used as a regularization parameter.'), ...
    arg({'ptype','Type'}, 'classification', {'classification','regression'}, 'Type of problem to solve.'), ...
    arg({'shape','Shape'}, [], [], 'Shape of the feature matrices. If given as [X,NaN] or [NaN,X], such that X is a divisor of the number of features F, the NaN is replaced by F/X.','shape','row'), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'), ...
    ... % misc
    arg({'vectorize_trials','VectorizeTrials'}, false, [], 'Vectorize trial matrices. Auto-determined if left unspecified.'), ...
    ... % parameter setup
    arg({'setupfcn','SetupFunction'}, [], [], 'Setup Function. Receives trials, targets, shape and the options structure, and returns the variables [X,y,B,pot,tau,G]. If empty, performs ridge regression or logistic regression depending on Type parameter.'), ...
    ... % inference engine parameters
    arg({'outerNiter','OuterIterations'},10,uint32([1 5 25 1000]),'Outer loop iterations.'),...
    arg({'outerMVM','OuterLanczosVectors'},50,uint32([1 10 100 1000]),'Outer Lanczos vectors. Number of Lanczos vectors to compute; more may yields a more complete posterior estimation.'),...
    arg({'innerMVM','InnerMVMSteps'},50,uint32([1 10 100 1000]),'Inner-loop MVM steps. Number of matrix-vector multiplications and/or conjugate gradient steps to perform in inner loop.'),...
    arg({'innerIt','InnerNewtonSteps'},15,uint32([1 5 50 1000]),'Inner-loop MVM steps. Number of Newton steps to perform in inner loop (if Truncated-Newton is used as the solver).'),...
    arg({'outerZinit','InitialUpperBound'},0.05,[],'Initial upper bound.','guru',true),...
    arg({'outerGainit','InitialVariationalParam'},1,[],'Initial variational parameter.','guru',true),...
    arg({'outerExact','ExactSolution'},false,[],'Exact solution. Whether the posterior covariance should be exactly computed (note: this may be prohibitively slow).'),...
    arg({'outerOutput','OuterVerbose','verbose'},false,[],'Verbose outer loop. Show diagonstic output in the outer loop of the algorithm.'),...
    arg({'innerOutput','InnerVerbose'},false,[],'Verbose inner loop. Show diagonstic output in the inner loop of the algorithm.'),...
    arg({'innerType','InferenceEngine'},'VariationalBayes',{'VariationalBayes','ExpectationPropagation'},'Inference engine. Inference engine to use in the inner loop; Variational Bayes is highly recommended as it works with all options, whereas Expectation Propagation comes with some restrictions (see glm-ie help).'),...
    arg({'innerVBpls','Solver'},'Conjugate Gradients',{'Quasi-Newton','Truncated Newton','Backtracking Conjugate Gradients','Conjugate Gradients','Split-Bregman'},'PLS solver. Penalized least-squares solver to use in the inner loop; Quasi-Newton (L-BFGS) is preferred but depends on Fortran code, (backtracking) Conjugate Gradients and truncated Newton are among the most useful fallbacks.'),...
    arg({'innerEPeta','InnerEPEta'},1,[],'Fractional EP parameter. Parameter used for fractional expectation-propagation (only works with Gauss and Laplace potential functions).','guru',true), ...
    arg({'doinspect','InspectMode'},false,[],'Inspection Mode. If enabled, the execution will break after the weights have been learned.','guru',true));

arg_toworkspace(args);

% set default setup function
if isempty(setupfcn)
    setupfcn = @default_setup; end

% get the correct feature matrix shape
if isempty(shape) %#ok<*NODEF>
    if ndims(trials) == 3
        shape = [size(trials,1) size(trials,2)];
        % ... also make sure that the trials are vectorized
        trials = double(reshape(trials,[],size(trials,3))');
        vectorize_trials = true;
    else
        shape = [size(trials,2) 1];
    end
elseif size(shape,1) == 1
    nf = size(trials,2);
    ni = isnan(shape);
    if any(ni)
        % if necessary, set NaN shape parameters appropriately
        shape(ni) = nf / shape(~ni);
    elseif nf ~= shape(1)*shape(2)
        % otherwise check for consistency
        error('shape parameter is inconsistent with feature space dimension.');
    end
end


% pre-process the data
if strcmp(ptype,'classification')
    classes = unique(targets);
    if length(classes) > 2
        % in the multi-class case we use the voter
        model = ml_trainvote(trials, targets, '1v1', @ml_trainglm, @ml_predictglm, varargin{:},'shape',shape,'vectorize_trials',vectorize_trials);
        return
    elseif length(classes) == 1
        error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
    else       
        % optionally scale the data
        sc_info = hlp_findscaling(trials,scaling);
        trials = hlp_applyscaling(trials,sc_info);
        % remap target labels to -1,+1
        targets(targets==classes(1)) = -1;
        targets(targets==classes(2)) = +1;
    end
elseif strcmp(ptype,'regression')
    classes = [];
    % scale the data
    sc_info = hlp_findscaling(trials,scaling);
    trials = hlp_applyscaling(trials,sc_info);
else
    error('Unrecognized problem type.');
end

% rewrite arguments
args.innerType = hlp_rewrite(args.innerType,'VariationalBayes','VB','ExpectationPropagation','EP');
args.innerVBpls = hlp_rewrite(args.innerVBpls,'Quasi-Newton','plsLBFGS','Truncated Newton','plsTN','Backtracking Conjugate Gradients','plsCGBT','Conjugate Gradients','plsCG','Split-Bregman','plsSB','Barzilai/Borwein','plsBB');

% call setup function
if nargout(setupfcn) == 5
    [X,y,B,pot,tau] = setupfcn(trials,targets,shape,args);
    G = [];
else
    [X,y,B,pot,tau,G] = setupfcn(trials,targets,shape,args);
end

% run inference
[uinf,ga,b,z,nlZ,Q,T] = hlp_diskcache('predictivemodels',@dli,X,y,lambdas,B,pot,tau,rmfield(args,{'trials','targets','lambdas','ptype','scaling','setupfcn','shape','doinspect'}),G); %#ok<ASGLU>

% add misc meta-data to the model
model.ptype = ptype;
model.classes = classes;
model.sc_info = sc_info;
model.shape = shape;
model.w = uinf;
model.vectorize = vectorize_trials;

% graphics
if doinspect
    disp('GLM inspection breakpoint; halted.');
    keyboard;
end


function [X,y,B,pot,tau,G] = default_setup(trials,targets,shape,opts)
[n,f] = size(trials);
if strcmp(opts.ptype,'regression')
    % ridge regression
    X = trials';
    y = targets;
    B = eye(f);
    pot = @potGauss;
    tau = ones(f,1);
    G = [];
else
    % logistic regression
    X = eye(f);
    y = zeros(f,1);
    B = trials';
    pot = @potLogistic;
    tau = targets;
    G = [];
end

