function [AR PE] = mvar_dalSCSA(varargin)
% Algorithm: Group Lasso (DAL/SCSA)
% 
% Description:
%
% This method infers a fully connected multivariate
% This program infers a sparsely connected multi-
% variate AR (VAR) model [1] assuming that sources 
% are independent and super Gaussian (hyperbolic 
% secant likelihood). An option is also provided to
% use a gaussian likelihood (squared loss function)
%
% VAR[p] coefficients inferred using Group Lasso 
% (L1,2 regularization) via the Sparsely Connected 
% Sources Analysis (SCSA) and Dual Augmented 
% Lagrangian (DAL) methods [1][2].
% The p VAR coefficients describing interactions 
% between two processes at time lags [t-1...t-p]
% are grouped together using an L2 norm which 
% penalizes large coefficients. An L1 penalty is
% then applied to the coefficient groups. This
% jointly prunes entire groups of coefficients
% by setting the entire group to zero. The result
% is a connectivity graph with sparse structure
% (most processes are strictly non-interacting)
% while regularizing (smoothing) the coefficient
% sequences describing surviving non-zero 
% interactions.
% These constraints allow us to uniquely solve 
% highly under-determined systems (e.g. many 
% more parameters than observations).
%
% If Y = D(y,p) is a p-lag delay embedding of multi-
% variate data vector y(t), X is a sparse block-
% diagonal matrix of lagged copies of the delay 
% embedded data, and A is an augmented matrix of
% VAR coefficients, then we may adopt the structural 
% model:
%
% Y = XA + e, for (super)gaussian noise e
%
% We then seek to solve the optimization problem
%
% [A_hat mu] = argmin_A{f(Y-XA) + L*sum(||A_i||_2)}
%
% Where f(z) specifies the unregularized portion of 
% the cost function (e.g. mimimizing variance of 
% residuals, e) and is given by:
% f(z) = ||z||_2^2  for gaussian innovations, and
% f(z) = -sum(log((1/pi)sech(z))) for supergaussian
% (sech) innovations.
% 
% where A_i contains the i^th group of AR parameters 
% e.g. the set of VAR weights {A(i,j)} describing auto-
% regression (conditional linear dependence) of X_i 
% onto X_j. L is the regularization parameter.
%
% This algorithm is a good choice if you have few 
% data samples and many channels/sources and/or a 
% high model order, or if the innovations are assumed
% nongaussian (as with ICA sources).
%
% Author Credits:
% 
% Adapted from s_test_hsgl.m from DAL 1.05 [3]
% Copyright(c) 2009-2011 Ryota Tomioka
%              2009      Stefan Haufe
%
% Reference:
%
% [1] S. Haufe, R. Tomioka, G. Nolte, K-R. Mueller, M. Kawanabe. IEEE Trans. Biomed. Eng. 57(8), pp. 1954-1963, 2010.
% [2] R. Tomioka, T. Suzuki, and M. Sugiyama. JMLR, 2011.
% [3] DAL Toolbox: http://www.ibis.t.u-tokyo.ac.jp/ryotat/dal/
%
% Dependencies: dalhsgl()
%
% ------------------------------------------------------------------------
% INPUTS:
%   data:       the data (nchs x npnts)
%   p:          the model order
%   lambda:     regularization coefficient
%   AR0:        initial solution (default: zeros)
%   varargin:   'name',value paired options to pass to dal.m
% OUTPUTS:
%   AR:         [nchs x nchs*p] VAR coefficient matrix
%   PE:         [nchs x nchs] estimated noise covariance matrix

% Modified by Tim Mullen, 2011 from s_test_hsgl.m from DAL 1.05
% - Modified for use with SIFT
% - Optimized significantly
% - Added support for multiple realization
% - Added residual covariance matrix estimation
%
% This program infers a sparsely connected multivariate AR
% model assuming that the sources are independent and super
% Gaussian. We use the hyperbolic secant likelihood.
%
% In our IEEE TBME paper, we also estimate mixing matrix
% through an EM algorithm. This is omitted here for the sake
% of simplicity. Thus, all state variables are directly measured.
%
% Reference:
% Modeling sparse connectivity between underlying brain sources
% for EEG/MEG. Stefan Haufe, Ryota Tomioka, Guido Nolte,
% Klaus-Robert Mueller, and Motoaki Kawanabe, IEEE
% Trans. Biomed. Eng. 57(8), pp. 1954-1963, 2010.
%
% Copyright(c) 2009-2011 Ryota Tomioka
%              2009      Stefan Haufe
% This software is distributed under the MIT license. See license.txt

persistent initAR;


g = arg_define(varargin, ...
                arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts].'), ...
                arg({'morder','ModelOrder','p'},10,[],'VAR Model order'), ...
                arg_nogui({'AR0','InitialState'},[],[],'Initial VAR coefficient matrix','shape','matrix','type','expression'), ...
                arg_sub({'dal_args','DAL_Options'},[],{...
                        arg({'lambda','RegularizationParam'},1,[0 Inf],'Regularization parameter'), ...
                        arg({'shrink_diagonal','ShrinkDiagonal'},true,[],'Apply sparsity constraint to diagonals (self-connections)'), ...
                        arg({'loss','LossFunction'},'HyperbolicSecant',{'HyperbolicSecant','Squared'},'Optimization loss function'), ...
                        arg({'verbosityLevel','Verbosity','verb'},0,{int32(0) int32(1) int32(2)},'Verbosity Level. 0 = basic, 2 = detailed.','type','int32'), ...
                        arg({'dalnvps','DAL_NVP_Args'},'{}','','Name-value pairs for DAL','type','expression','shape','row')}, ...
                    'Options for DAL algorithm.') ...
                );

% 'maxiter',10,'tol',1e-3
            
arg_toworkspace(g);
p = morder;

[nchs npnts ntr] = size(data);

% initial solution
if ~isempty(AR0)
    initAR = AR0;
elseif isempty(initAR)
    initAR = zeros(p*nchs^2,1);
end
if size(initAR,1) ~= p*nchs^2
    % dimensions have changed, reset state
    initAR = zeros(p*nchs^2,1);
end

X = [];
Y = [];
for itr = 1:ntr
    Xi = [];
    for jj = 1:p
        Xi = cat(3, Xi, squeeze(data(:, p+1-jj:end-jj, itr))');
    end
    X = cat(1, X, reshape(permute(Xi, [1 3 2]), npnts-p, nchs*p));
    Y = cat(1,Y,data(:, p+1:end,itr)');
end

Y = vec(Y);
    

% ne    = ntr*(npnts-p);            % number of block equations of size m
% Mp    = nchs*p;
% 
% X = zeros(ne,Mp);
% Y = zeros(ne,nchs);
% 
% for itr=1:ntr
%     for j=1:p
%         % Build matrix of predictors (right hand side of regression model)
%         X((npnts-p)*(itr-1) + 1 : (npnts-p)*itr , nchs*(j-1)+1 : nchs*j) = ...
%             squeeze(data(:, p-j+1:npnts-j, itr))';
%     end
%     % Build observation matrix (left hand side of regression model)
%     Y((npnts-p)*(itr-1) + 1 : (npnts-p)*itr, 1 : nchs) = squeeze(data(:, p+1:npnts, itr))';
% end

X = blkdiageye(sparse(X),nchs);



%% Indices for the diagonal elements
diagIdx = vec(repmat((1:p*(nchs+1):p*nchs^2), p, 1) + repmat((0:(p-1))', 1, nchs))';

%% Indices for the off-diagonal elements
offDiagIdx = setdiff_bc(1:p*nchs^2, diagIdx);

%% Design corresponding to the diagonal elements (self connection)
B = X(:, diagIdx);
if dal_args.shrink_diagonal
    blks = [p*ones(1,nchs*(nchs-1)),p*nchs];
    Bu = [];
else
    Bu = B;
    B = [];
    blks = p*ones(1,nchs*(nchs-1));
end

%% Design corresponding to the off-diagonal elements (connection to others)
X = X(:, offDiagIdx);

%% Target
Y = Y(:);
% Y = vec(data(:, p+1:end)');

% Group penalize AR coefficients with different time-lags
% between sources (nchs*(nchs-1) groups of size p).
% Self connections are put into one group (size p*nchs).
if dal_args.shrink_diagonal
    bias = []; %zeros(size(initAR));
    
    switch lower(dal_args.loss)
        case 'squared'
            initAR=dalsqgl(initAR, [X B], Y, dal_args.lambda,  ...
                           'blks', blks,'verbosityLevel',dal_args.verbosityLevel, ...
                           dal_args.dalnvps{:});
        case 'hyperbolicsecant'
            initAR = dalhsgl(initAR,bias, [X B],[], ...
                             Y, dal_args.lambda, 'blks', blks, ...
                             'verbosityLevel',dal_args.verbosityLevel, ...
                             dal_args.dalnvps{:});  % 'aa',zeros(size(X,1),1),
    end
else
    [initAR(offDiagIdx) initAR(diagIdx)] = dalhsgl(initAR(offDiagIdx),initAR(diagIdx), [X B],Bu, ...
                                                   Y, dal_args.lambda, 'blks', blks, ...
                                                   'verbosityLevel',dal_args.verbosityLevel,dal_args.dalnvps{:});
end

H2 = zeros(p,nchs,nchs);
H2(offDiagIdx) = vec(initAR(1:(p*nchs*(nchs-1))));          % non-diagonal elements
H2(diagIdx) = vec(initAR(((p*nchs*(nchs-1))+1):end));       % diagonal elements

AR = permute(full(H2), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
AR = reshape(AR,[nchs nchs*p]);

%% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(data,AR,zeros(1,nchs));
    res = res(:,:);
    PE = cov(res',1);
end

% if nargout>2
%     argsout.initAR = initAR;
% end

function v = vec(x)
v = x(:);


function C = blkdiageye(X,k)
% construct blockdiagonal matrix with k copies of X on diagonal
% Equivalent to C = kron(eye(k),X) but much faster

ss = repmat('X,',1,k); ss(end)=[];
eval(sprintf('C = blkdiag(%s);',ss));

