
function Rinv = est_calcInvCovMat(AR,C,invCov,mode,MAXITER)
%
% For any M-variate VAR[p] process fit to data X = [X1, X2, ... Xn], 
% this function returns the [(Mp)^2 x (Mp)^2] covariance (or inverse covariance) 
% matrix of the process. This is a block-toeplitz, hermitian matrix with 
% format as described in [2]. 
% 
% Briefly, let M=nvars, and Rhat be comprised of M x M submatrices 
%
%           Rhat(u,v) = [R11(u,v) ... R1M
%                           ...        ...
%                        RM1(u,v) ... RMM(u,v)]
%
% where Rij(u,v) = cov(Xi(t-u),Xj(t-v)) for i,j=1:M and u,v = 1:p. Then
% Rhat is the (cross-)covariance matrix of the VAR[p] process up to lag p
% with inverse covariance matrix Rinv = inverse(Rhat).
%
% Inputs:
%
%   AR:     [nvars x nvar*morder] matrix of vector autoregressive
%           coefficients
%   C:      nvars x nvars noise covariance matrix
%   invCov: {def: 1}
%           true: return inverse of processes covariance matrix
%           fase: return process covariance matrix
%   mode:   {def: 3} Determines the method used to compute the covariance matrix
%           1: uses the standard approach from [2]
%           2: same as 1, but using sparse matrices.
%           3: uses the iterative doubling algorithm of Anderson and Moore
%              as proposed in [1]. For large matrices, generally much faster
%              and less memory intensive.
%   MAXITER:{def: 100}
%           Number of iterations for doubling algorithm (mode=3)
%
% Outputs:
%
%   Rinv:   [(nvars*p)^2 x (nvars*p)^2] (inverse) covariance matrix of the
%           VAR[p] process.
%
%           
% References:
%
% [1] Barone, (1987) J. of Time series Analysis vol. 8 no. 2
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
% Springer. Sec 2.1.4.
% 
% See Also: est_calcInvCovMatFourier(), est_calcInvCovMatPDC
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD
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



if nargin<4
    mode = 3;
end
if nargin<3
    invCov = 1;    % if 0, don't invert covariance matrix
end
if nargin<4
    MAXITER = 100;
end

% make sure the noise covariance matrix is valid
C = covfixer(C);

switch mode
    
    case 1    
        
        % Standard covariance estimation. See [2] sec. 2.1.4
        
        M = size(C,1);          % nchans
        p = size(AR,2)/M;       % model order
        Mp=M*p;

        A = zeros(Mp);
        A(1:M,1:Mp) = AR;
        A(M+1:end,1:end-M) = eye(Mp-M);
        Eu = zeros(Mp);
        Eu(1:M,1:M) = C;

        Rinv = (eye(Mp^2)-kron(A,A))\Eu(:);

%         Rinv = reshape(Rinv,Mp,Mp)\eye(Mp);
        
        Rinv = inverse(reshape(Rinv,Mp,Mp));

    case 2    
        
        % same as 1, but sparse approach
        
        M = size(C,1);          % nchans
        p = size(AR,2)/M;       % model order
        Mp=M*p;

        A = spalloc(Mp,Mp,M*Mp+M);
        A(1:M,1:Mp) = AR;
        A(M+1:end,1:end-M) = eye(Mp-M);

        Eu = spalloc(Mp,Mp,M^2);
        Eu(1:M,1:M) = single(C);

        Rinv = (speye(Mp^2)-kron(A,A))\Eu(:);

        % Rinv = symmlq(speye(Mp^2)-kron(A,A),Eu(:),1e-6,100);

        Rinv = full(reshape(Rinv,Mp,Mp)\speye(Mp));
        
    case 3
        
        % Anderson-Moore/Barone Doubling Algorithm [1]

        M = size(C,1);          % nchans
        p = size(AR,2)/M;       % model order
        Mp=M*p;
        A = zeros(Mp);
        A(1:M,1:Mp) = AR;
        A(M+1:end,1:end-M) = eye(Mp-M);

        K = [eye(M); zeros(M*(p-1),M)];
        N = {K*C*K'};
        N{2} = zeros(size(N{1}));
        
        for t=1:MAXITER
            N{2} = A*N{1}*A'+N{1};
            A = A^2;
            if norm(N{2}(:)-N{1}(:))<1e-15
                break
            end
            N{1}=N{2};
        end
        
        if t==MAXITER
            warning('SIFT:est_calcInvCovMat', 'Anderson-Moore did not converge. Results may be innacurate.\n');
        end
        
        % Rinv = L*N{2}*L';

        % ensure that the covariance matrix is a valid covariance matrix
        N{2} = covfixer(N{2});
            
        if invCov
            Rinv = inverse(N{2});
        else
            Rinv = N{2};
        end
    
end
