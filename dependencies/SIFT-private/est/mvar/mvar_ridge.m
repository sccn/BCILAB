function [AR PE lambdaOpt] = mvar_ridge(varargin)
% Algorithm: Ridge Regression
%
% Description:
%
% This method infers a fully connected multivariate
% autoregressive (VAR) model using L2 regularization,
% otherwise known as "ridge-regression", "Tikhonov 
% regularization" or "minimum-L2-norm" [2]
%
% VAR[p] coefficients inferred using ridge regression
% (L2 regularization). We assume most VAR coefficients
% (interactions) are small, and apply an L2 penalty
% for large coefficients. A consequence is that VAR
% coefficients are never exactly zero; thus statistical
% thresholding is an important follow-up step.
% This constraint allows us to solve under-determined
% systems (e.g. more parameters than observations).
%
% If Y = D(y,p) is a p-lag delay embedding of multi-
% variate data vector y(t), X is a block toeplitz 
% matrix of lagged copies of the delay embedded data,
% and A is an augmented matrix of VAR coefficients,
% then we may adopt the structural model:
%
% Y = XA + e,  for gaussian noise e
%
% We then seek to solve the optimization problem
%
% A_hat = argmin_A{0.5||Y-XA||_2^2 + L*||A||_2}
%
% where L is the regularization parameter. 
% L can be determined automatically using Generalized 
% Cross Validation (GCV) [3].
%
% This algorithm is a good choice if you have few 
% data samples and a fairly large number of sources
% and/or a high model order. Note that, unlike sparse
% methods, ridge regression yields a fully connected
% graph (all VAR coefficients are non-zero). Follow up
% with statistical thresholding.
%
% Author Credits:
% 
% Implemented by Tim Mullen.
% The ridge regression implementation (ridgeGCV.m) was
% contributed by Alejandro Ojeda (SCCN/INC).
%
% Dependencies: ridge_gcv()
%
% ------------------------------------------------------------------------
% INPUTS:
%   data:       the data (nchs x npnts)
%   p:          the model order
%   lambda:     regularization coefficient
% OUTPUTS:
%   AR:         [nchs x nchs*p] VAR coefficient matrix
%   PE:         [nchs x nchs] estimated noise covariance matrix
%
% See Also: est_fitMVAR(), mvar_dalSCSA(), est_fitMVARKalman()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual.
%     Available at: http://www.sccn.ucsd.edu/wiki/SIFT
% [2] Hoerl, A.E.; R.W. Kennard (1970). Ridge regression: Biased estimation 
%     for nonorthogonal problems. Technometrics 42(1):80?86. JSTOR 1271436.
% [3] Golub G, Heath M, Wahba G. (1979) Generalized Cross-validation as a  
%     method for choosing a good ridge parameter. Technometrics. 
%     doi: 10.1080/00401706.1979.10489751
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

verb = arg_extract(varargin,{'verb','Verbosity'},[],0);

if hlp_isToolboxInstalled('Parallel Computing Toolbox')
    pardef   = 'on';
    try
    [tmp parprofs] = hlp_microcache('sift_domain',@defaultParallelConfig);
    catch err
        pardef = 'off';
        parprofs = {'local'};
    end
else
    pardef = 'off';
    parprofs = {'local'};
end

g = arg_define(varargin, ...
                arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts x ntr].'), ...
                arg({'morder','ModelOrder','p'},10,[],'VAR Model order'), ...
                arg({'normcols','NormCols'},'norm',{'none','norm','zscore'},'Normalize data. Columns of data matrix X and target vector Y are normalized using the chosen method.'), ...
                arg_subtoggle({'splitVars','DecoupleVars','SplitVars'},'on', ...
                { ...
                    arg_subtoggle({'runPll','RunInParallel'},'off', ...
                    { ...
                    arg({'profile','ProfileName'},parprofs{1},parprofs,'Profile name'), ...
                    arg({'numWorkers','NumWorkers'},2,[1 Inf],'Number of workers') ...
                    },'Solve blocks in parallel. Requires Parallel Computing Toolbox.'), ...
                },'Decouple variables. Solve coefficients for each variable separately. This can improve computation speed when there are a large number of channels/variables'), ...
                arg_sub({'ridge_args','RegressionOptions'},[],@ridge_gcv,'Ridge regression options.','suppress','verb'), ...
                arg({'verb','Verbosity'},verb,{int32(0) int32(1) int32(2)},'Verbose output','type','int32','mapper',@(x)int32(x)) ...
                );

if strcmp(pardef,'off') && g.splitVars.arg_selection && g.splitVars.runPll.arg_selection
    fprintf('Parallel Computing Toolbox not installed. Cannot use parallel option.\n');
    g.splitVars.runPll.arg_selection = false;
end

arg_toworkspace(g);

[nchs npnts ntr] = size(data);
p = morder;

% assemble the predictor (design) matrix and target vector
[X Y] = hlp_mkVarPredMatrix(g.data,p);

% design matrix block size (one block per variable)
blkrows   = npnts-p;
blkcols   = p*nchs;
blksz = [blkrows*ntr, blkcols];

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
X = blkdiageye(sparse(X),nchs);


if g.splitVars.arg_selection
    
    AR = zeros(size(X,2),1);
    lambdaOpt  = zeros(1,nchs);
    
    % solve for each variable separately
    if g.splitVars.runPll.arg_selection
        if g.verb
            fprintf('Running parallel job (this may take a while)...\n'); 
        end
        ridge_args = g.ridge_args;
        
        % make sure pool is open
        wasOpen = hlp_pll_openPool(g.splitVars.runPll.numWorkers, ...
                                   g.splitVars.runPll.profile);
        
        % slice the variables
        for k=1:nchs
            rowidx = (k-1)*blksz(1)+1:k*blksz(1);
            colidx = (k-1)*blksz(2)+1:k*blksz(2);
            Xk{k}  = X(rowidx,colidx);
        end
        Yk = reshape(Y,blksz(1),nchs);
        verb = g.verb;
        % run first iteration to obtain svd_state
        [ARk{1}, lambdaOpt(1), svd_state] = ridge_gcv('Y',Yk(:,1),'A',Xk{1}, ...
                                                     ridge_args,'verb',verb);
        parfor k=2:nchs
            [ARk{k}, lambdaOpt(k)] = ridge_gcv('Y',Yk(:,k),'A',Xk{k}, ...
                                                     ridge_args,    ...
                                                    'svd_state',svd_state,'verb',verb);
        end
        % reassemble results
        for k=1:nchs
            AR((k-1)*blksz(2)+1:k*blksz(2)) = ARk{k};
        end
        % cleanup
        if ~wasOpen
            matlabpool('close'); end
    else
        svd_state = [];
        for k=1:nchs
            rowidx = (k-1)*blksz(1)+1:k*blksz(1);
            colidx = (k-1)*blksz(2)+1:k*blksz(2);
            [AR(colidx), lambdaOpt(k), svd_state] ...
                = ridge_gcv('Y',Y(rowidx),'A',X(rowidx,colidx),g.ridge_args,   ...
                            'blksz',[],'svd_state',svd_state,'verb',g.verb);
        end
    end
    
else
    [AR lambdaOpt] = ridge_gcv('Y',Y,'A',X,g.ridge_args,'blksz',blksz,'verb',g.verb);
end

AR = reshape(AR,[p,nchs,nchs]);
AR = permute(full(AR), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
AR = reshape(AR,[nchs nchs*p]);

%% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(data,AR,zeros(1,nchs));
    res = res(:,:);
    PE  = cov(res',1);
end
