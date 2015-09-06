function [X Y] = hlp_mkVarPredMatrix(data,p)
% Return the predictor (design) matrix, X, and response matrix, Y, of the 
% order-1 (VAR[1]) equivalent model for a zero-mean VAR[p] model.
%
% Suppose we have data, x, of dimension [M x N x T] where 
% M = nchs  = number of channels
% N = npnts = number of data samples per trial
% T = ntr   = number of trials
% 
% For simplicity, assume T = ntr = 1 (the case T>1 is handled automatically
									  % in this function, and discussed later).
%
% The zero-mean, order-p, vector autoregression equation
%
% x(:,t) = A1*x(:,t-1) +...+ Ap*x(:,t-p) + noise(:,t)     [Eq. 1]
%
% can be rewritten as an order-1 VAR model:
%
% Y = X*A + noise     [Eq. 2]
% 
% where Y is the 'response' (target) matrix, X is the 'predictor' (design)
% matrix and A contains the unknown VAR coefficients to be estimated:
% 
% X     = [X1 ... Xp]
% Xk    = [x(:, p+1-k) ... x(:, N-k)]'
% Y     = X0
% A     = [A1 ... Ap]'
%
% X is of size [(N-p) x M*p]  ( [T*(N-p) x M*p] for T>1 )
% Y is of size [(N-p) x M]    ( [T*(N-p) x M]   for T>1 )
% A is of size [M*p x M]
%
% Many solvers can exploit known structure in the model (e.g. Block SBL, 
% Group Lasso, etc). One such structure is the natural grouping of the p
% coefficients {Aij} that describe influences from channel j to channel i
% across all lags. For convenience when exploiting structure, we want to 
% reorder the columns of X (rows of A) such that coefficients within a group 
% are stored in consecutive rows of A and the corresponding predictors are 
% stored in consecutive columns of X.
%
% With this modification, we have the following matrix design:
%
% X     = [X1 ... XM]
% Xi    = [x(i, p+1-1)' ... x(i, N-p)']
%
% X is a block-Toeplitz matrix consisting of M blocks [X1 ... XM] where
% each block Xi is a Toeplitz matrix of size [npnts-p x p] whose p columns 
% contain the time-delayed data for channel i. Specifically, the kth column 
% of Xi is the (p-k+1)th-order delay embedding of the data for channel i.
% That is, Xi(:,k) is the data vector for channel i delayed by p-k+1
% samples (smaller k => greater delay).
%
% ----- Multi-trial extension (T>1) -----
%
% If T>1, we simply treat each trial as an independent observation of the 
% predictor and response matrices for the VAR[1] model specified in [Eq 2]. 
% As such, the T-trial augmented regression model can be written as a 
% system of T VAR[1] regression models of the form in [Eq. 2]:
%
% Y1  =  X1*A + e1
% Y2  =  X2*A + e2
%    ...
% YT  =  XT*A + eT
% 
% which is equivalent to the augmented VAR[1] model:
%
% | Y1 |      | X1 |            | e1 |
% | Y2 |  =   | X2 | * [A]  +   | e2 |
% | ...|      | ...|            | ...|
% | YT |      | XT |            | eT |
% -------------------------------------
%   Y     =     X    *  A   +    noise
%
% For each trial, we form predictor/response matrices as described above. 
% We then vertically stack these matrices to form the appropriate
% multi-trial X and Y (note this does not change the dimensions of A). 
%
% ----- Note on parallelizing the solution to Y = XA + noise -----
% 
% Note that the solution to the regression problem
% 
% Y = X*A
% 
% is given by 
% 
% pinv(X)*Y = A 
% 
% where pinv() is some pseudoinverse operator. We can see that the ith col.
% of A is obtained by projecting the ith column of Y onto the inverse of X.
%
% pinv(X)*Y(:,i) = A(:,i)
%
% As such, we can solve for each column of A independently using the 
% corresponding column of Y yeilding M independent models of the form:
% 
% Y(:,i) = X*A(:,i) + noise
% 
% which can solved in M parallel operations.
%
% -------------------------------------------------------------------------
% Inputs:
%   data:   [nchs x npnts x ntr] data matrix
%   p   :   model order
% Ouputs:
%   X   :   [ntr*(npnts-p) x nchs*p] predictor matrix
%   Y   :   [ntr*(npnts-p) x nchs]   target matrix

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


[nchs npnts ntr] = size(data);

% get the dimensions of the delay-embedding blocks
blkrows   = npnts-p;
blkcols   = p*nchs;

% Assemble X, Y
X = zeros(blkrows,blkcols,ntr);
Y = zeros(blkrows,nchs,ntr);
for itr = 1:ntr
    % initialize delay-embedding blocks
    Xi = zeros(blkrows,nchs,p);
    for d = 1:p % extract data at each delay lag...
        Xi(:,:,d) = data(:, p+1-d:end-d, itr)';
    end
    % permute to [npnts x p x nchs] so each page is a 
    % delay-embedding block for a given channel...
    Xi = permute(Xi, [1 3 2]); 
    % ... and concatenate blocks horizontally to form 2D design mat
    % for this trial X(:,:,itr) = [X1 X2 ... XM]
    X(:,:,itr) = reshape(Xi, blkrows, blkcols);
    
    % extract 0-lag data matrix for this trial...
    Y(:,:,itr) = data(:, p+1:end,itr)';
end

if ntr>1
    % reshape X and Y to stack trials vertically and form final 2D matrix
    X = permute(X,[1 3 2]);
    Y = permute(Y,[1 3 2]);
    X = reshape(X,blkrows*ntr,blkcols);
    Y = reshape(Y,blkrows*ntr,nchs);
end

