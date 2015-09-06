function [AR PE State] = mvar_BSBL_L1(varargin)
% Algorithm: BSBL L1
%
% Description:
%
% BSBL Algorithm by Zhilin Zhang
%
% References and Code:
%
% Dependencies: BSBL_L1_noise()
%
% ------------------------------------------------------------------------
% INPUTS:;
%   data:       the data (nchs x npnts)
%   p:          the model order
%   lambda:     regularization parameter
%   rho:        the augmented Lagrangian parameter
%   alpha:      the over-relaxation parameter (typical 
%               values for alpha are between 1.0 and 1.8)
%   initState:        initial solution (default: zeros)
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
%       Michael Jordan, Editor in Chief, 3(1):1?122, 2011.
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

persistent initAR;

g = arg_define(varargin, ...
                arg_norep({'data','Data'},mandatory,[],'Data Matrix. Dimensions are [nchs x npnts].'), ...
                arg({'morder','ModelOrder'},10,[],'VAR Model order'), ...
                arg({'normcols','NormCols'},'zscore',{'none','norm','zscore'},'Normalize columns of dictionary'), ...
                arg_subtoggle({'warmStart','WarmStart'},[], ...
                {...
                    arg({'initState','InitialState'},[],[],'Initial BSBL state object. If empty, defaults will be used.') ...
                },'Warm start. The previously estimated state will be used as a starting estimate for successive BSBL operations. Alternately, you may provide an non-empty initial state structure via the ''InitialState'' argument.'), ...
                arg({'groupDiags','GroupAutoConnections','GroupDiags'},false,[],'Group auto-connections. All auto-connections for all channels will be penalized jointly as a separate group. Note, this can slow down model estimation as the design matrix will not longer be diagonal block-toeplitz so we cannot (yet) exploit block-redundancy in the design matrix'), ...
                arg_sub({'bsbl_args','BSBL_Options'},[],{ ...
                       arg({'LearnLambda','LambdaLearningRule'},'LowNoise',{'off','LowNoise','HighNoise'},{'Lambda learning rule',sprintf(['\n' ...
                                                                                                              '''LowNoise'' :  Use the lambda learning rule for very LOW SNR cases (SNR<10dB). Uses lambda=std(y)*1e-2 or user-input value as initial value.\n\n' ...
                                                                                                              '''HighNoise'':  Use the lambda learning rule for medium noisy cases (SNR>10dB). Uses lambda=std(y)*1e-2 or user-input value as initial value.\n\n' ...
                                                                                                              '''off'':        Do not use the lambda learning rule. Uses lambda=1e-14 or user-input value as initial value' ...
                                                                                                              ])}) ...
                       arg({'PRUNE_GAMMA','GammaPruningThreshold'},[],[0 Inf],{'Gamma Pruning threshold',sprintf(['\n' ...
                                                                                   'Threshold for prunning small hyperparameters gamma_i. Blocks of parameters are pruned when their ''power'' (gamma_i) is small.\n'    ...
                                                                                   'In noisy cases, you can set PRUNE_GAMMA = 1e-3 or 1e-4. \n'   ...
                                                                                   'In strong noisy cases (e.g. SNR <= 6 dB), set PRUNE_GAMMA = 1e-2 for better performance.' ...
                                                                                 ])}), ...
                       arg({'LAMBDA','InitialLambda','Lambda'},[],[],'Initial regularization value (lambda). If empty, the initial value is estimated from the data.','type','denserealdouble'), ...                                                        
                       arg({'EPSILON','StoppingTolerance'},1e-8,[0 Inf],'Stopping tolerance. The algorithm terminates when the L2 norm of the change in parameter estimates is small than this value.'), ...
                       arg({'LEARNTYPE','LearnCorrelation'},true,[],'Exploit correlation structure.'), ...
                       arg({'MAX_ITERS','MaxIterations'},10,[1 Inf],'Maximum iterations'), ...
                       arg({'PRINT','VerboseOutput'},false,[],'Verbosity') ...
                   },'Additional options for BSBL algorithm') ...
            );
                
% arg_toworkspace(g);
  
% modify some arguments
switch g.bsbl_args.LearnLambda
    case 'off',         g.bsbl_args.LearnLambda = 0;
    case 'LowNoise',    g.bsbl_args.LearnLambda = 1;
    case 'HighNoise',   g.bsbl_args.LearnLambda = 2;
end

[nchs npnts ntr] = size(g.data);
p = g.morder;
blkrows   = npnts-p;

% initialize state
if g.warmStart.arg_selection && ~isempty(g.warmStart.initState)
    initAR = g.warmStart.initState;
elseif ~g.warmStart.arg_selection
    % reset initAR
    initAR = zeros(p*nchs^2,1);
end
if size(initAR,1) ~= p*nchs^2
    % dimensions have changed, reset state
    fprintf('mvar_BSBL_L1: model dimensions changed -- resetting state\n');
    initAR = zeros(p*nchs^2,1);
end

% -------------------------------------------------------------------------
% Build the predictor (design) matrix and target vectors
% -------------------------------------------------------------------------

% assemble the predictor (design) matrix and target vector
[X Y] = hlp_mkVarPredMatrix(g.data,p);

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
% % Apply the BSBL method for sparse VAR estimation:
% -------------------------------------------------------------------------
% group penalize AR coefficients with different time-lags
% between sources (nchs*(nchs-1) groups of size p). 
g.bsbl_args.InitState = initAR;
args = hlp_struct2varargin(g.bsbl_args,'suppress',{'arg_direct','LearnLambda'});
[initAR] = BSBL_L1_noise(X,Y,blks,g.bsbl_args.LearnLambda,args{:});


% assemble coefficient matrices
AR = zeros(p,nchs,nchs);
if g.groupDiags
    AR(offDiagIdx) = vec(initAR(1:(p*nchs*(nchs-1))));          % non-diagonal elements
    AR(diagIdx) = vec(initAR(((p*nchs*(nchs-1))+1):end));       % diagonal elements
else
    AR(:) = initAR;
end

AR = permute(full(AR), [3 2 1]);   % the recovered (nchs x nchs x p) connectivity matrix
AR = reshape(AR,[nchs nchs*p]);

%% estimate noise covariance matrix
if nargout>1
    res = est_mvarResiduals(g.data,AR,zeros(1,nchs));
    res = res(:,:);
    PE = cov(res',1);
end

if nargout>2
    State = initAR;
end

function v = vec(x)
v = x(:);

