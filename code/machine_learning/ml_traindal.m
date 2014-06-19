function model = ml_traindal(varargin)
% Learn a linear probabilistic model via the Dual-Augmented Lagrangian method.
% Model = ml_traindal(Trials, Targets, Lambda, Options...)
%
% The Dual-Augmented Lagrangian method [1] is an efficient and very robust approach to learning
% regularized linear classifiers or regressors, particularly for "noisy" biosignals. The
% regularization is very effective, so that a large number of features (e.g. every channel and time
% point) can be supplied for learning. Assumptions about how features are correlated or independent
% w.r.t. each other can (and should) be incorporated, by specifying the appropriate type of
% regularizer. If features are known to be all uncorrelated (e.g. derived from indepenent
% components), 'l1' is the appropriate regularizer. If features are correlated only within blocks,
% 'glr'/'glc' (group lasso by rows or columns) is the appropriate regularizer, e.g. time points of
% concatenated independent components. If features are all correlated, but it is understood that
% there when arranged in a feature matrix, the correlation structure would be low-rank, then 'ds' is
% the appropriate regularizer. The 'ds' mode is ideally suited when features are a matrix of
% time-points by channels, where both time-points are mutually correlated and channels are, too, or
% when the features are covariance matrices, etc.
%
% To inform the classifier of the block size for 'glr'/'glc' or the matrix shape for 'ds', the
% trials should either be supplied as a 3d matrix (i.e. feature matrices instead of feature
% vectors), or supplied in the regular fashion (2d matrix of feature vectors), but with the intended
% feature matrix shape specified in the 'shape' option. A few methods with different
% performance/accuracy/datasize tradeoffs are supplied. An important consideration when using DAL is
% that the data must be appropriately normalized for the method to be most effective, that is, it
% must be normalized across features and/or groups in 'l1' and 'glc'/'glr' modes, and normalized
% across both the horizontal and vertical axes of the feature matrix, and/or across blocks of a
% block-diagonal feature matrix, in the 'ds' mode.
%
% BCI paradigms which make extensive use of this classifier, according to [2], are provided in the
% paradigms/para_dal* functions. Among the methods provided in the toolbox, DAL is likely the best
% applicable method if the data is linearly separable (albeit not necessarily the easiest to use).
%
% In:
%   Trials       : training data, as in ml_train
%                  in addition, it may be specified as UxVxN 3d matrix,
%                  with UxV-formatted feature matrices per trial (N trials), or
%                  as {{U1xV1,U2xV2,...}, {U1xV1,U2xV2,...}
%
%   Targets      : target variable, as in ml_train
%
%   Lambda       : sequence of regularization parameters to evaluate; (default: 2.^(10:-1:-15))
%
%   Options  : optional name-value parameters to control the training details:
%              'loss': loss function to be used,
%                        'squared' for regression
%                        'logistic' for classification (default)
%
%              'regularizer': type of regularization to use:
%                             'l1': l1-norm on the features, gives sparse results
%                             'glr': grouped l1 norm ("group lasso"), gives blockwise sparse results
%                                   (groups the rows of the feature matrices)
%                             'glc': grouped l1 norm ("group lasso"), gives blockwise sparse results
%                                   (groups the columns of the feature matrices)
%                             'ds': dual spectral norm, gives low-rank results (default)
%                             'en': elastic net norm, employs a combination of l1 and l2 regularization
%
%              'shape': if trials is a NxF 2d matrix of vectorized matrices,
%                           this is the dimensions of the matrices (default: Fx1)
%                       if trials is specified as an UxVxN 3d matrix, shape defaults to
%                           [U,V] and trials are vectorized into the regular [N, U*V]
%                       if shape is specified with one of the values being NaN,
%                           that value is set so that prod(shape)==F==U*V
%
%              misc parameters:
%              'scaling': pre-scaling of the data (see hlp_findscaling for options) (default: 'none')
%
%              'nfolds' : Cross-validation folds. The cross-validation is used to determine the best 
%                         regularization parameter value (default: 5)
%
%              'foldmargin' : Margin (in trials) between folds. This is the number of trials omitted 
%                             between training and test sets. (default: 5)
%
%              'cvmetric' : metric to use for parameter optimization; can be any of those supported by
%                           ml_calcloss (default: '' = auto-determine)
%
%              'bias': whether to include a bias term (default: 1)
%
%              'quiet': whether to suppress diagnostic outputs (default: 1)
%
%              'solver': solver to be used:
%                         'cg'   : Newton method with preconditioned conjugate gradient descent (default)
%                         'qn'   : Quasi-Newton method (slower, but uses less flimsy code...)
%
% Out:
%   Models   : a predictive model
%
% Examples:
%   % assuming a 3d feature array of size UxVxT, and a label vector of size Tx1
%   % the features should be appropriately normalized (see [2] for examples)
%
%   % learn a DAL model for a given regularization parameter using the logistic loss (for classification)
%   model = ml_traindal(trials,targets,0.1)
%
%   % as before, but use the squared loss, for linear regression
%   model = ml_traindal(trials,targets,0.1,'loss','squared')
%
%   % like before, but this time use the 'l1' (LASSO) regularizer, assuming sparse features
%   model = ml_traindal(trials,targets,0.1,'regularizer','l1')
%
%   % like before, but this time use the group LASSO regularizer imposing group sparsity on the rows
%   % of the feature matrix (columns is 'glc')
%   model = ml_traindal(trials,targets,0.1,'regularizer','glr')
%
%   % like before but use the (default) dual-spectral regularizer, which learns low-rank weights
%   model = ml_traindal(trials,targets,0.1,'regularizer','ds')
%
%   % if the individual trials are not matrix-shaped but vectorized, pass in the shape manually
%   model = ml_traindal(trials,targets,0.1,'shape',[U,V])
%
%   % use a different solver (here: conjugate gradient, which is potentially more efficient)
%   model = ml_traindal(trials,targets,0.1,'solver','cg')
%
%   % learn a DAL model using a parameter search
%   model = utl_searchmodel({trials,targets},'args',{{'dal', search(2.^(10:-0.5:-6))}})
%   
%   % as before, but use a different loss
%   model = utl_searchmodel({trials,targets},'args',{{'dal', search(2.^(10:-0.5:-6)), 'loss', 'glc'}})
%
% See also:
%   ml_predictdal, dal
%
% References:
%  [1] Ryota Tomioka & Masashi Sugiyama, "Dual Augmented Lagrangian Method for Efficient Sparse Reconstruction",
%      IEEE Signal Proccesing Letters, 16 (12) pp. 1067-1070, 2009.
%  [2] Ryota Tomioka and Klaus-Robert Mueller, "A regularized discriminative framework for EEG analysis with application to brain-computer interface",
%      Neuroimage, 49 (1) pp. 415-432, 2010.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-25

arg_define([0 3],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'lambdas','Lambdas'}, 2.^(4:-0.25:-3), [0 2^-5 2^15 Inf], 'Regularization parameters. Controls the sparsity/simplicity of the result. Typically, this is an interval to scan, such as 2.^(10:-1:-15).'), ...
    arg({'loss','LossFunction'}, 'logistic', {'logistic','squared'}, 'Loss function . The logistic loss is suited for classification problems, whereas the squared loss is suited for regression problems.'), ...
    arg({'regularizer','Regularizer'}, 'dual-spectral', {'lasso','grouplasso-rows','grouplasso-columns','dual-spectral'}, 'Type of regulariation to use. Lasso (l1) gives sparse results (e.g., on ic-spectral decompositions), the grouped l1 norms give blockwise sparse results (rows / columns of the feature matrices), e.g. for spatially or spectrally decomposed data, and the dual-spectral norm gives low-rank results (raw eeg)'), ...
    arg({'shape','Shape'}, [], [], 'Shape of the feature matrices. If given as [X,NaN] or [NaN,X], such that X is a divisor of the number of features F, the NaN is replaced by F/X.','shape','row'), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'), ...
    arg({'nfolds','NumFolds'},5,[0 Inf],'Cross-validation folds. The cross-validation is used to determine the best regularization parameter.'),...
    arg({'foldmargin','FoldMargin'},5,[0 0 10 Inf],'Margin between folds. This is the number of trials omitted between training and test set.'), ...    
    arg({'theta','Theta'}, 0.5, [0 1], 'Elastic net blending. Blend parameter for the elastic net.'), ...    
    arg({'cvmetric','ParameterMetric'},'',{'','kld','nll','mcr','mae','mse','max','rms','smse','sign','bias','medse','auc','cond_entropy','cross_entropy','f_measure'},'Metric for Parameter Optimization. By default auto-determined; can be any of the ml_calcloss-supported metrics. In particular, auc is a good idea if the classification task is between highly imbalanced classes.'), ...
    arg({'verbose','Verbose'},false,[],'Show diagnostic output.'), ...
    arg({'solver','Solver'},'conjugate gradient',{'conjugate gradient','quasi-Newton'},'Solution method. These differ in robustness, speed and memory requirements.'),...
    arg({'votingScheme','VotingScheme'},'1v1',{'1v1','1vR'},'Voting scheme. If multi-class classification is used, this determine how binary classifiers are arranged to solve the multi-class problem. 1v1 gets slow for large numbers of classes (as all pairs are tested), but can be more accurate than 1vR.'), ...
    arg_deprecated({'doinspect','InspectMode'},false,[],'Inspection Mode. Not used any more -- use a breakpoint instead.'));


classes = unique(targets);
if length(classes) > 2 && strcmp(loss,'logistic')
    % in the multi-class case we use the voter
    model = ml_trainvote(trials, targets, votingScheme, @ml_traindal, @ml_predictdal, varargin{:});
elseif length(classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
else
    if strcmp(cvmetric,'mcr')
        cvmetric = ''; end
    
    % get the correct feature matrix shape
    vectorize_trials = false;
    if isempty(shape) %#ok<*NODEF>
        if ndims(trials) == 3
            shape = [size(trials,1) size(trials,2)];
            % ... also make sure that the trials are vectorized
            trials = double(reshape(trials,[],size(trials,3))');
            vectorize_trials = true;
        else
            shape = [size(trials,2) 1];
            if any(strcmp(regularizer,{'dual-spectral','ds','grouplasso-rows','glr','grouplasso-columns','glc'}))
                warn_once('BCILAB:DAL:ill_advised_usage','You are using the DAL method with a regularizer that makes sense only on group-structured features, which are however not specified. Falling back to lasso. Consider using logreg in the LARS variant instead.');
                regularizer = 'lasso';
            end
        end
    elseif size(shape,1) == 1
        nf = size(trials,2);
        ni = isnan(shape);
        if any(ni)
            % if necessary, set NaN shape parameters appropriately
            shape(ni) = nf / shape(~ni);
        elseif length(shape) >= 3
            % tensor-shaped features are grouped over 2d tensor pages
            shape = repmat([shape(1),shape(2)],prod(shape(3:end)),1);
        elseif nf ~= shape(1)*shape(2)
            % otherwise check for consistency
            error('shape parameter is inconsistent with feature space dimension.');
        end
    end
    
    if isempty(lambdas)
        lambdas = 2.^(10:-0.25:-5); end
    
    % optionally scale the data
    sc_info = hlp_findscaling(trials,scaling);
    trials = hlp_applyscaling(trials,sc_info);
    
    % rewrite the bias, regularizer & solver to the format expected by DAL
    regularizer = hlp_rewrite(regularizer,'lasso','l1','grouplasso-rows','glr','grouplasso-columns','glc','dual-spectral','ds','elastic-net','en');
    solver = hlp_rewrite(solver,'conjugate gradient','cg','quasi-Newton','qn','Newton with Cholesky decomposition','nt','Newton with memory saving','ntsv','subspace trust-region','fminunc');
    
    % remap target labels to -1,+1
    if strcmp(loss,'logistic')
        targets(targets==classes(1)) = -1;
        targets(targets==classes(2)) = +1;
    end
    
    % possibly the data needs to be transposed
    dotranspose = strcmp(regularizer,'glr');
    if dotranspose
        shape = shape([2 1]);
        trials = double(reshape(trials',shape(2),shape(1),[]));
        ntrials = zeros(shape(1),shape(2),size(trials,3));
        for t=1:size(trials,3)
            ntrials(:,:,t) = trials(:,:,t)'; end
        trials = double(reshape(ntrials,[],size(ntrials,3))');
    end
    
    % lambdas need to be sorted in descending order for the warm-starting to work...
    lambdas = sort(lambdas,'descend');

    % determine the correct learning function to use, according to loss & regularizer...    
    switch loss
        case 'logistic'
            switch regularizer
                case 'ds'
                    learner = @dallrds; % dual-spectral logistic regression
                case {'glc','glr'}
                    learner = @dallrgl; % group-regularized logistic regression
                case 'l1'
                    learner = @dallrl1; % sparse logistic regression
                case 'en'
                    learner = @(ww,bias,A,yy,lambda,varargin) dallren(ww,bias,A,yy,lambda,theta,varargin{:}); % elastic-net logistic regression
                otherwise
                    error('Unsupported regularizer.');
            end
        case 'squared'
            switch regularizer
                case 'ds'
                    learner = @dalsqds; % dual-spectral regularized regression
                case {'glc','glr'}
                    learner = @dalsqgl; % group LASSO
                case 'l1'
                    learner = @dalsql1; % LASSO
                case 'en'
                    learner = @(ww,A,bb,lambda,varargin) dallren(ww,A,bb,lambda,theta,varargin{:}); % elastic net
                otherwise
                    error('Unsupported regularizer.');
            end
        otherwise
            error('Unsupported loss function.');
    end
    
    % learn an ensemble of models across the given lambda's, on all the data (i.e. the regularization path)
    if verbose
        disp('Running DAL...'); end

    if length(lambdas) > 1        
        % cross-validate to score the lambda's
        foldid = 1+floor((0:length(targets)-1)/length(targets)*nfolds); 
        reptargets = repmat(targets,1,length(lambdas));
        predictions = zeros(length(targets),length(lambdas));
        % for each fold...
        for i = 1:nfolds
            if verbose
                disp(['Fitting fold # ' num2str(i) ' of ' num2str(nfolds)]); end
            
            % determine training and test set indices
            which = foldid==i;
            trainids = ~which;
            whichpos = find(which);
            for j=1:foldmargin
                trainids(max(1,whichpos-j)) = false;
                trainids(min(length(which),whichpos+j)) = false;
            end
            
            % learn an ensemble of models...
            subensemble = hlp_diskcache('predictivemodels',@learn_ensemble,learner,lambdas,shape,trials(trainids,:),targets(trainids),solver,loss,verbose,regularizer);
            
            % obtain test-set predictions for each model...
            testset = [trials(which,:) ones(length(whichpos),1)];
            for m=1:length(subensemble)
                curmodel = subensemble{m};
                w = full([curmodel.w(:); curmodel.b]);
                predictions(which,m) = (testset*w)';
            end
            if strcmp(loss,'logistic')
                predictions(which,:) = 2*(1 ./ (1 + exp(-predictions(which,:))))-1; end

            % evaluate the loss
            if isempty(cvmetric)
                if strcmp(loss,'logistic')
                    loss_mean(i,:) = mean(reptargets(which,:) ~= sign(predictions(which,:)));
                else
                    loss_mean(i,:) = mean((reptargets(which,:) - predictions(which,:)).^2);
                end
            else
                for r=1:length(lambdas)
                    loss_mean(i,r) = ml_calcloss(cvmetric,reptargets(which,r),predictions(which,r)); end
            end
        end
        loss_mean = mean(loss_mean,1);

        % legacy losses output
        if strcmp(loss,'logistic')
            losses = reptargets ~= sign(predictions);
        else
            losses = (reptargets - predictions).^2;
        end
        
        % if there are several minima, choose largest lambda of the smallest cvm
        lambda_min = max(lambdas(loss_mean <= min(loss_mean)));
    else
        lambda_min = lambdas;
        losses = NaN;
    end
    
    % pick the model at the minimum...
    ensemble = hlp_diskcache('predictivemodels',@learn_ensemble,learner,lambdas,shape,trials,targets,solver,loss,verbose,regularizer);    
    model = ensemble{find(lambdas == lambda_min,1)};
    model.transpose = dotranspose;
    model.classes = classes;
    model.sc_info = sc_info;
    model.shape = shape;
    model.loss = loss;
    model.ensemble = ensemble;
    model.losses = losses;
    model.vectorize = vectorize_trials;
end



% learn the regularization path using the DAL method...
function ensemble = learn_ensemble(learner,lambdas,shape,trials,targets,solver,loss,verbose,regularizer)
% learn_ensemble_version<1.0.0>
disp('learning ensemble...');

% derive the design matrix A & label vector y from the trials...
if size(shape,1) == 1
    mask = any(trials);
    if mean(mask) < 0.75
        [m,n] = size(trials);
        % check if block-diagonal
        tmask = double(trials(:,mask));
        fA = @(x) tmask * x(mask);
        fAt = @(x) spreadvec(tmask'*x,mask,n);
        A = {fA,fAt,m,n};
    elseif nnz(trials)/numel(trials) < 0.25
        % check if reasonably sparse
        A = sparse(trials);
    else
        A = double(trials);
    end
else
    A = double(trials);
end
y = double(targets);

% now learn the models
ensemble = cell(1,length(lambdas));
curmodel = struct('w',{zeros(sum(shape(:,1).*shape(:,2)),1)},'b',{0});
for k =1:length(lambdas)
    lam = lambdas(k);
    fprintf('  scanning lambda = %f...',lam);
    % scale up lambda by nTrials since the data term is not normalized by 1/nTrials
    lam = lam * length(targets);
    % learn an updated model
    t0 = tic;
    if strcmp(loss,'logistic')
        if any(strcmp(regularizer,{'glc','glr'}))
            [curmodel.w,curmodel.b] = learner(reshape(curmodel.w(:),shape),curmodel.b,A,y,lam,'display',verbose,'solver',solver);
        else
            [curmodel.w,curmodel.b] = learner(curmodel.w(:),curmodel.b,A,y,lam,'display',verbose,'solver',solver,'blks',shape);
        end
    else
        curmodel.w = learner(curmodel.w(:),A,y,lam,'display',verbose,'solver',solver,'blks',shape);
    end
    duration = toc(t0);
    try
        % calculate the final rank...
        ix = 0;
        modelrank = 0;
        for s=1:size(shape,1)
            ival = shape(s,1)*shape(s,2);
            modelrank = modelrank + rank(reshape(curmodel.w(ix+1:ix+ival),shape(s,:)));
            ix = ix+ival;
        end
        % display diagnostics
        fprintf(' model rank = %i; t = %.1fs\n',modelrank,duration);
    catch
        fprintf('\n');
    end
    % store
    ensemble{k} = curmodel;
end


% spread a sparse vector out according to an index set
function y = spreadvec(x,idx,n)
y = zeros(n,1);
y(idx) = x;
