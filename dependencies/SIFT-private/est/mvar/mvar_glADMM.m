function [AR PE] = mvar_glADMM(varargin)

% Algorithm: Group Lasso (ADMM)
%
% Description:
%
% This method infers a sparsely connected multi-
% variate autoregressive (VAR) model under a 
% gaussian noise assumption.
%
% VAR[p] coefficients inferred using Group Lasso 
% (L1,2 regularization) via the Alternating 
% Direction Method of Multipliers [1]. 
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
% Y = XA + e,  for gaussian noise e
%
% We then seek to solve the optimization problem
%
% A_hat = argmin_A{0.5||Y-XA||_2^2 + L*sum(||A_i||_2)}
%
% where A_i contains the i^th group of AR parameters 
% e.g. the set of VAR weights {A(i,j)} describing auto-
% regression (conditional linear dependence) of X_i 
% onto X_j. L is the regularization parameter.
%
% This algorithm is a good choice if you have few 
% data samples and many channels/sources and/or a 
% high model order.
%
% Author Credits:
%
% This implementation is based on Matlab ADMM
% examples from Stephen Boyd's website [2]
%
% References and Code:
%
% [1] Boyd, Parikh, Chu, Pelato and Eckstein et al. Foundations and Trends in Machine Learning 3(1):1-122,2011
% [2] http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
% Dependencies: admm_gl()
%
% ------------------------------------------------------------------------
% INPUTS:
%   data:       the data (nchs x npnts)
%   p:          the model order
%   lambda:     regularization parameter
%   rho:        the augmented Lagrangian parameter
%   alpha:      the over-relaxation parameter (typical 
%               values for alpha are between 1.0 and 1.8)
%   AR0:        initial solution (default: zeros)
%   verb:       true/false - verbosity
% OUTPUTS:
%   AR:         [nchs x nchs*p] VAR coefficient matrix
%   PE:         [nchs x nchs] estimated noise covariance matrix
%
% See Also: est_fitMVAR(), mvar_dalSCSA(), mvar_ridge(), mvar_vieiramorf(),
%           est_fitMVARKalman()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
% [2] S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein, "Distributed 
%       Optimization and Statistical Learning via the Alternating Direction 
%       Method of Multipliers". Foundations and Trends in Machine Learning, 
%       Michael Jordan, Editor in Chief, 3(1):1-122, 2011.
%       http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
% [3] http://www.stanford.edu/~boyd/papers/admm/group_lasso/group_lasso_example.html
%
% Author: Tim Mullen 2012, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

persistent initAR gpusel;

% set some arg defaults
% GPUDeviceNames = {};
% if hlp_microcache('mvar_glAdmm',@hlp_isToolboxInstalled,'Parallel Computing Toolbox')
%      [GPUDeviceNames, GPUDeviceIdx, GPUDeviceCC] = hlp_microcache('mvar_glADMM',@hlp_getGPUDeviceNames,1.3);
%      % append compute capability to device name
%      %GPUDeviceNames = cellfun(@(dn,cc) sprintf('%s_CC%0.3g',dn,cc),GPUDeviceNames,GPUDeviceCC);
%      if isempty(GPUDeviceNames)
%          GPUDeviceNames = {'none'};
%      else
%         GPUDeviceNames = [{'default'} GPUDeviceNames];
%      end
% else
%     GPUDeviceNames = {'default'};
% end
GPUDeviceNames = {'default'};
% parse inputs
g = arg_define(varargin, ...
                arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts].'), ...
                arg({'morder','ModelOrder','p'},10,[],'VAR Model order'), ...
                arg_subtoggle({'warmStart','WarmStart'},[], ...
                {...
                    arg({'initState','InitialState'},[],[],'Initial ADMM state. This is a structure with fields ''z_init'' and ''u_init'', which represent, respectively, the initial state vector and dual vector and are both of dimension [morder*(nchs^2) x 1]. If empty, the first state is initialized to zeros.') ...
                },'Warm start. The previously estimated state will be used as a starting estimate for successive operations. Alternately, you may provide an non-empty initial state structure via the ''InitialState'' argument.'), ...
                arg({'normcols','NormCols'},'none',{'none','norm','zscore'},'Normalize columns of dictionary'), ...
                arg({'groupDiags','GroupAutoConnections','GroupDiags'},false,[],'Group auto-connections. All auto-connections for all channels will be penalized jointly as a separate group. Note, this can slow down model estimation as the design matrix will not longer be diagonal block-toeplitz so we cannot (yet) exploit block-redundancy in the design matrix'), ...
                arg_norep({'AR0','InitialState'},[],[],'DEPRECATED. Initial VAR coefficient matrix','shape','matrix','type','expression'), ...
                arg_subtoggle({'gpu','GPU'},'off',...
                {
                    arg({'device','Device'},GPUDeviceNames{1},GPUDeviceNames,'Which GPU device to use') ...
                    arg({'library','Library'},'CUDA_ADMM',{'CUDA_ADMM'},'Which library to use'), ...
                },'Use GPU. This uses the CUDA_ADMM package by Oleg Konigs.'), ...    
                arg_sub({'admm_args','ADMM_Options'},[],@admm_gl,'Options for ADMM algorithm') ...
                );

% arg_toworkspace(g);
if g.gpu.arg_selection && strcmp(g.gpu.device,'none')
    error('No GPU device available with sufficient compute capability');
end
        
if g.gpu.arg_selection && strcmp(g.gpu.library,'CUDA_ADMM') && ~exist('glADMM_CUDA','file')
    error('CUDA_ADMM package not installed. Please download and install from https://github.com/OlegKonings/GLwithmex');
end

% determine whether A can be compressed to a single block
can_optimize_X = ~g.groupDiags && (~g.gpu.arg_selection && strcmp(g.admm_args.x_update.arg_selection,'direct'));

[nchs, npnts, ntr] = size(g.data);
p = g.morder;
blkrows   = npnts-p;
blkcols   = p*nchs;

if isempty(initAR)
    initAR.z = [];
    initAR.u = [];
end

% initialize state
if g.warmStart.arg_selection && ~isempty(g.warmStart.initState)
    initAR = g.warmStart.initState;
elseif ~g.warmStart.arg_selection
    % reset initAR
    initAR.z = zeros(p*nchs^2,1);
    initAR.u = zeros(p*nchs^2,1);
end
if size(initAR.z,1) ~= p*nchs^2
    % dimensions have changed, reset state
    fprintf('mvar_glADMM: model dimensions changed -- resetting state\n');
    initAR.z = zeros(p*nchs^2,1);
    initAR.u = zeros(p*nchs^2,1);
end

% -------------------------------------------------------------------------
% Build the predictor (design) matrix and target vectors
% -------------------------------------------------------------------------

% assemble the predictor (design) matrix and target vector
[X, Y] = hlp_mkVarPredMatrix(g.data,p);

% normalization
switch lower(g.normcols)
    case 'zscore'
        % standardize columns of X and Y (unit variance)
        X = zscore(X,1,1);
        Y = zscore(Y,1,1);
    case 'norm'
        % normalize columns of X and Y (unit norm)
        X = bsxfun(@rdivide,X,sqrt(sum(X.^2)));
        Y = bsxfun(@rdivide,Y,sqrt(sum(Y.^2)));
end

% It is convenient to express A as a column vector so we must...
% vectorize the target matrix...
Y = hlp_vec(Y);

% ...and expand predictor matrix
if ~can_optimize_X
    X = blkdiageye(sparse(X),nchs);
end

% If auto-prediction coefficients (connectivity matrix 'diagonals') are to
% be smoothed separately, we place them in a special block at the end of
% the coefficient vector. To acheive this, we also need to move the
% associated predictors for these coefficients to the rightmost columns
% of the predictor matrix
if g.groupDiags
    % Indices for the diagonal elements (self connection)
    diagIdx = vec(repmat((1:p*(nchs+1):p*nchs^2), p, 1) + repmat((0:(p-1))', 1, nchs))';
    % Indices for the off-diagonal elements (connection to others)
    offDiagIdx = setdiff_bc(1:p*nchs^2, diagIdx);

    X = [X(:, offDiagIdx) X(:, diagIdx)];
    blks = [p*ones(1,nchs*(nchs-1)),p*nchs];
else
    blks = p*ones(1,nchs^2);
end
    
% -------------------------------------------------------------------------
% Apply the ADMM method for group lasso estimation:
% -------------------------------------------------------------------------
% group penalize AR coefficients with different time-lags
% between sources (nchs*(nchs-1) groups of size p). 
if g.gpu.arg_selection
    % use GPU implementation
    
    % if a new GPU device has been selected, activate it
    if ~strcmp(g.gpu.device,'default') && ~strcmpi(gpusel,g.gpu.device)
        gpuDevice(GPUDeviceIdx(strcmp(g.gpu.device,GPUDeviceNames)));
        gpusel = g.gpu.device; % cache the device name
    end
    switch g.gpu.library
        case 'CUDA_ADMM'
%             glADMM_CUDA(AA,b,partition,u,z,rho,alpha,lambda,MAX_ITER,ABSTOL,RELTOL)
            % call the CUDA library
            % first make sure dimensions of X,Y,u,z are multiples of 16
            [~, n] = size(X);
            if (n > 1), colpad = 16*(floor(n/16)+logical(mod(n,16)))-n;
            else        colpad = 0; end
            [initAR.u,initAR.z]= glADMM_CUDA(...
                                     multpad(single(full(X')),16),      ...
                                     multpad(single(Y),16),             ...
                                     int32([blks colpad]'),             ...
                                     multpad(single(initAR.u),16),      ...
                                     multpad(single(initAR.z),16),      ...
                                     single(g.admm_args.rho),           ...
                                     single(g.admm_args.alpha),         ...
                                     single(g.admm_args.lambda),        ...
                                     int32(g.admm_args.max_iter),       ...
                                     single(g.admm_args.abstol),        ...
                                     single(g.admm_args.reltol));
            % trim excess zeros (if added during padding)
            if numel(initAR.u)>p*nchs^2, initAR.u(p*nchs^2+1:end) = []; end
            if numel(initAR.z)>p*nchs^2, initAR.z(p*nchs^2+1:end) = []; end
    end
else
    % use CPU implementation
    [initAR.z, initAR.u] = admm_gl('A',X,  ...
                           'y',Y, ...
                           'blks',blks, ...
                           g.admm_args, ...
                           'z_init',initAR.z, ...
                           'u_init',initAR.u, ...
                           'designMatrixBlockSize',fastif(g.groupDiags,[],[blkrows*ntr blkcols]), ...
                           'arg_direct',true);
end
% assemble coefficient matrices
AR = zeros(p,nchs,nchs);
if g.groupDiags
    AR(offDiagIdx) = vec(initAR.z(1:(p*nchs*(nchs-1))));          % non-diagonal elements
    AR(diagIdx)    = vec(initAR.z(((p*nchs*(nchs-1))+1):end));    % diagonal elements
else
    AR(:) = initAR.z;
end

AR = permute(full(AR), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
AR = reshape(AR,[nchs nchs*p]);

%% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(g.data,AR,zeros(1,nchs));
    res = res(:,:);
    PE = cov(res',1);
end

function v = vec(x)
v = x(:);


function X = multpad(X,base)
% pad dims of X with zeros such that dimensions are a multiple of 'base'
% only works for vectors and 2-D matrices

[m, n] = size(X);
if m==1 && n==1
    warning('scalar input provided to multpad.\nCannot determine expansion dimensions');
    return;
end
% determine the number of elements to pad rows and columns to enforce
% multiple of base
if (m > 1), rowpad = base*(floor(m/base)+logical(mod(m,base)))-m;
else        rowpad = 0; end
if (n > 1), colpad = base*(floor(n/base)+logical(mod(n,base)))-n;
else        colpad = 0; end

% pad the data
if rowpad, X = cat(1,X,zeros(rowpad,size(X,2),class(X))); end
if colpad, X = cat(2,X,zeros(size(X,1),colpad,class(X))); end
