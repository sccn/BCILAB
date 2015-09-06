function [C R scale dof] = est_qrNoiseCovMat(data,morder,mcor)
% Compute an estimate of the noise covariance matrix of an AR model
% using the QR decomposition approach of [1]. This program based on 
% code from ARFIT [2][3].
%
% References:
%
% [1] A. Neumaier and T. Schneider. ACM TOMS, 27(1):27-57, 2001.
% [2] T. Schneider and A. Neumaier. ACM TOMS, 27(1):58-65, 2001.
% [3] http://www.gps.caltech.edu/~tapio/arfit/
%
% Author: Tim Mullen 2013, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu


% n:   number of time steps (per realization)
% m:   number of variables (dimension of state vectors) 
% ntr: number of realizations (trials)
[n,m,ntr] = size(data);  
ne        = ntr*(n-morder);         % number of block equations of size m
np        = m*morder+mcor;          % maximum number of parameter vectors of length m

% compute QR factorization for model of order morder
[R, scale]   = arqr(data, morder, mcor);

% get lower triangle R22 according to 
%
%   | R11  R12 |
% R=|          |
%   | 0    R22 |
%  
R22   = R(np+1:np+m, np+1:np+m);

% compute covariance matrix
dof   = ne-np;                % number of block degrees of freedom
C     = R22'*R22./dof;        % bias-corrected estimate of covariance matrix

