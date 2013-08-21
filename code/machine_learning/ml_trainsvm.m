function model = ml_trainsvm(varargin)
% Learn a predictive model by Support Vector Machines.
% Model = ml_trainsvm(Trials, Targets, Cost, Options...)
%
% Support vector machines are nowadays some of the most frequently used and versatile linear
% classifiers [1]. They can handle a large number of feature efficiently and are usually extremely
% robust (if regularized well) or fast (when not regularized). The main differences between the
% provided variants are computational, and allow to select the most efficient model for the data (in
% terms of #trials/#features). All variants except for 'native' are implemented using the LIBLINEAR
% package [2], the 'native' variant servers as fallback (in case LIBLINEAR is not available for some
% platform). Support vector machines do not naturally produce good probability estimates (although
% some approach exists, like in the LDA case).
%
% in typical uses of the toolbox, SVMs are not as frequently used as logreg or lda, primarily due to
% the need for proper regularization, which is often prohibitively costly. For quick tests, the
% method can be used without regularization, however. When extremely many feature or trials are
% used, SVMs are likely more useful than in standard settings.
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%
%   Cost         : regularization parameter, reasonable range: 2.^(-5:2:15), greater is stronger
%                  default: 1
%
%   Options  : optional name-value pairs; possible names are:
%              'variant': one of several SVM variants
%                       'dual'    : L2-regularized L2-loss support vector classification (dual), usually preferred (default)
%                       'primal'  : L2-regularized L2-loss support vector classification (primal), can be faster than dual
%                       'crammer' : Crammer and Singer approach to multi-class support vector classification,
%                                   as alternative to the voted classification employed in the other variants
%                       'sparse'  : L1-regularized L2-loss support vector classification, gives sparse results but
%                                   will likely not be better unless the solution is truly sparse.
%                       'l1loss'  : L2-regularized L1-loss support vector classification (dual), rarely needed
%                       'native'  : native MATLAB implementation (using CVX); equivalent to 'dual', but slower
%
%              'kernel': one of several kernel types:
%                         * 'linear':   Linear
%                         * 'rbf':      Gaussian / radial basis functions (default)
%                         * 'laplace':  Laplacian
%                         * 'poly':		Polynomial
%                         * 'cauchy':	Cauchy
%
%              'gamma': scaling parameter of the kernel (for regularization), default: 0.3
%                        reasonable range: 2.^(-16:2:4)
%
%              'degree': degree of the polynomial kernel, if used (default: 3)
%
%              'eps'    : desired accuracy (default: [], meaning LIBLINEAR's default)
%
%              'bias'   : include a bias variable (default: 1)
%
%              'scaling': pre-scaling of the data (see hlp_findscaling for options) (default: 'std') 
%
%              'verbose'  : verbosity level (default: 0)
%
% Out:
%   Model   : a predictive model; 
%             classes indicates the class labels which the model predicts
%             sc_info is the scaling info
%
% Examples:
%   % learn a quick and dirty SVM classifier (without parameter search)
%   model = ml_trainsvm(trials,labels)
%
%
%   % learn an SVM model by searching over the cost parameter (note: quite slow)
%   model = utl_searchmodel({trials,labels},'args',{{'svm',search(2.^(-5:2:15))}})
%
%   % as before, but also search over the kernel scale parameter (note: really slow)
%   model = utl_searchmodel({trials,labels},'args',{{'svm',search(2.^(-5:2:15)),'gamma',search(2.^(-16:2:4))}})
%
%   % as before, but use a linear kernel (no need to search over gamma in this case)
%   model = utl_searchmodel({trials,labels},'args',{{'svm',search(2.^(-5:2:15)),'kernel','linear'}})
%
%   % as before, but use primal optimization; this is faster under certain circumstances
%   model = utl_searchmodel({trials,labels},'args',{{'svm',search(2.^(-5:2:15)),'variant','primal','kernel','linear'}})
%
%   % as before, but use the native MATLAB implementation; this may be used when the binary code 
%   % does not run for some reason (note: should better be run on a cluster, as this is extremely slow)
%   model = utl_searchmodel({trials,labels},'args',{{'svm',search(2.^(-5:2:15)),'variant','native','kernel','linear'}})
%
%
% See also:
%   ml_predictsvm
%
% References:
%   [1] P. S. Bradley and O. L. Mangasarian. "Massive data discrimination via linear support vector machines."
%       Optimization Methods and Software, 13:1-10, 2000. 
%   [2] R.-E. Fan, K.-W. Chang, C.-J. Hsieh, X.-R. Wang, and C.-J. Lin. "LIBLINEAR: A library for large linear classification" 
%       Journal of Machine Learning Research 9(2008), 1871-1874.
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-04

arg_define([0 3],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'cost','Cost'}, search(2.^(-5:2:15)), [], 'Regularization parameter. Reasonable range: 2.^(-5:2:15), greater is stronger.'), ...
    arg({'variant','Variant'}, 'dual', {'dual','primal','crammer','sparse','l1loss','native'}, 'Variant to use. Dual and primal are l2-regularized l2-loss support vector classification, dual usually preferred but primal can be faster. Crammer is a special approach to multi-class problems. Sparse assumes that the result should be sparse, l1loss uses a rarely used loss function, native is a native MATLAB implementation, if the binary fails to work.'), ...
    arg({'kernel','Kernel'}, 'rbf', {'linear','rbf','laplace','poly','cauchy'}, 'Kernel type. Linear, or Non-linear kernel types: Radial Basis Functions (general-purpose), Laplace (sparse), Polynomial (rarely preferred), and Cauchy (slightly experimental).','cat','Core Parameters'), ...
    arg({'gammap','KernelScale','gamma'}, search(2.^(-16:2:4)), [], 'Scaling of the kernel functions. Should match the size of structures in the data. A reasonable range is 2.^(-16:2:4).','cat','Core Parameters'), ...
    arg({'polydegree','PolyDegree','degree'}, uint32(3), [], 'Degree of the polynomial kernel, if chosen.','cat','Core Parameters'), ...    
    arg({'epsi','Epsilon','eps'}, [], [], 'Tolerated solution accuracy. If unspecified, the default of the library (LIBLINEAR) will be taken.'), ...
    arg({'bias','Bias'}, true, [], 'Include a bias term.'), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'), ...
    arg({'verbose','Verbose'}, false, [], 'Whether to show diagnostic output.'));

if is_search(cost)
    cost = 1; end
if is_search(gammap)
    gammap = 0.3; end


% identify and remap the classes (to 0..k)
classes = unique(targets); %#ok<*NODEF>
if strcmp(variant,'native')
    % MATLAB version
    if length(classes) > 2
        % in this case we use the voter...
        model = ml_trainvote(double(trials),targets,'1v1',@ml_trainsvm,@ml_predictsvm,varargin{:});
        return;
    elseif length(classes) == 1
        error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
    end
    % scale the data
    sc_info = hlp_findscaling(trials,scaling);
    trials = hlp_applyscaling(trials,sc_info);
    % kernelize the data
    basis = trials;
    trials = utl_kernelize(trials,basis,kernel,gammap,polydegree);    
    % remap targets
    targets(targets == classes(1)) = -1;
    targets(targets == classes(2)) = +1;    
    % solve
    [n,f] = size(trials); %#ok<NASGU,ASGLU>
    cvx_begin
        variables w(f) b xi(n)
        minimize 1/2*sum(w.*w) + cost*sum(xi)
        subject to
            targets.*(trials*w + b) >= 1 - xi;
            xi >= 0;
    cvx_end
    model = struct('w',w,'b',b);
else
    % scale the data
    sc_info = hlp_findscaling(trials,scaling);
    trials = hlp_applyscaling(trials,sc_info);        
    % kernelize the data
    basis = trials;
    trials = utl_kernelize(trials,basis,kernel,gammap,polydegree);
    % remap targets
    targ = targets;
    for c=1:length(classes)
        targets(targ==classes(c)) = c-1; end
    
    % check whether we are using the fallback variant if LIBLINEAR
    fallback = ~isempty(strfind(which('lltrain'),'pre7.3'));
    if fallback && strcmp(variant,'sparse')
        error('Sparse classifiers are not supported in LIBLINEAR prior to MATLAB 7.3 (due to a breaking change in the sparse array MEX API).'); end    
    
    % rewrite some args
    variant = hlp_rewrite(variant,'dual',1,'primal',2,'crammer',4,'l1loss',3,'sparse',5);
    bias = hlp_rewrite(bias,true,'-B 1',false,'');
    quiet = hlp_rewrite(~verbose && ~fallback,true,'-q',false,'');
    
    % build the arguments
    if ~isempty(epsi)
        args = sprintf('-s %d -c %f -e %f %s %s',variant,cost,epsi,quiet,bias);
    else
        args = sprintf('-s %d -c %f %s %s',variant,cost,quiet,bias);
    end
    % run the command
    model = hlp_diskcache('predictivemodels',@llwtrain,ones(size(targets)),targets,sparse(double(trials)),args);
end
% finalize the model
model.sc_info = sc_info;
model.classes = classes;
model.variant = variant;
model.basis = basis;
model.kernel = kernel;
model.gammap = gammap;
model.degree = polydegree;