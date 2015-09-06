function [res]=est_mvarResiduals(X,AR,mu,downsampleFactor,instant)
%
%  Computes the M-variate residuals for variables k=1...M
%
%    res(:,k) = X(:,k+p) - mu - A1*X(:,k+p-1) - ... - Ap*X(:,k)
%
%  of an VAR(p) model with AR=[A1 ... Ap]
%
%   Inputs:
%
%       X:      [M x npoints x ntrials] matrix of sensor data
%       AR:     VAR[p] model coefficients in format AR=[A1 ... Ap]. Can be
%               a cell array of length (npnts) where AR{i} are the VAR
%               coefficients for the ith time point or window
%       mu:     [1 x nchannels] vector of model intercepts (if no intercepts
%               estimated, set mu = [] or mu = zeros(1,nchannels) and process mean
%               is assumed to be zero
%       downsampleFactor:   This is the downsample factor used in
%                           est_fitMVARKalman() or, equivalently,
%                           ceil(srate*winstep). Note that, since we are
%                           computing prediction errors, if the
%                           coefficients are time-varying, this function
%                           will return the d-step prediction error where
%                           d=downsampleFactor (the previous coefficient
%                           matrix -- which could correspond to many
%                           time-points in the past -- is used to predict the
%                           current timepoint). To obtain "instantaneous"
%                           prediction errors, set the instant option to true.
%       instant:            If true, compute "instantaneous" prediction errors
%                           (use the AR(t) to predict X(t)), else
%                           compute N-step prediction error (use AR(t-d) to
%                           predict X(t), where d = downsampleFactor)
%   Output:
%
%       res:    [M x npoints-p x ntrials] matrix of residuals (1-step ahead
%               prediction error). Note that
%
% See Also: pop_est_fitMVAR()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual. Chapter 4
%     Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD.
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

if ~iscell(AR)
    AR = {AR};
end

nchs  = size(X,1);                    % number of variables
n     = size(X,2);                    % number of observations
ntr   = size(X,3);                    % number of realizations (trials)
p     = size(AR{1},2)/nchs;           % order of model
tvar  = length(AR)>1;                 % boolean determining whether VAR is time-varying
if ~tvar
    nres  = n-p;                    % number of residuals
else
    nres  = length(AR);
end

if nargin<3 || isempty(mu)
    mu = zeros(nchs,1);
else
    mu     = mu(:);                     % force mu to be column vector
end
if nargin<4
    downsampleFactor = 1;
end
if nargin<5
    instant = true;
end

res = nan*ones(nchs,nres,ntr);

if ~tvar
    mu     = mu*ones(1,nres);          % construct matrix
    
    % Get time series of residuals
    l = 1:nres;
    % vectorized loop l=1,...,nres
    for tr=1:ntr
        res(:,l,tr) = X(:,l+p,tr) - mu;
        for k=1:p
            res(:,l,tr) = res(:,l,tr) - AR{1}(:, (k-1)*nchs+1:k*nchs)*X(:,l-k+p,tr);
        end
    end
    
    % reference code
%     tmp = zeros(nchs,nres,ntr);
%     for k=1:p
%         tmp = tmp + AR{1}(:, (k-1)*nchs+1:k*nchs)*X(:,l-k+p,tr);
%     end
% 
%     % code vectorization tests (not ready)
%     ARx = reshape(AR{1},nchs,nchs,p);
%     M = reshape(permute(ARx(:,:,end:-1:1),[1 3 2]),[],nchs) * X(:,:,tr);
%     idx = bsxfun(@plus,bsxfun(@plus,[1:nchs]',[0:nres-1]*(nchs*p)),permute([0:p-1]*(nchs*p+nchs),[1 3 2]));
%     N = reshape(M(idx(:)),nchs,p,[]);
%     tmp2 = reshape(sum(N,2),nchs,[]);
    
else
    if ~instant
        % add initial AR coefficients (this delays coefficient matrix by one
        % sample, which produces d-step prediction error)
        AR = [zeros(size(AR{1})) AR];
    end
    
    for tr=1:ntr
        idx = 1;
        for t=max(2,downsampleFactor):downsampleFactor:n
            if t<=p
                Y = [X(:,t-1:-1:1,tr) zeros(nchs,p-t+1)];
            else
                Y = X(:,t-1:-1:t-p,tr);
            end
            
            z = AR{ceil((t-1)/downsampleFactor)}.';
            
            H = blkdiageye(sparse(Y(:)'),nchs);
            res(:,idx,tr) = X(:,t)-mu-H*z(:);
            idx = idx+1;
        end
    end
end

%   % Center residuals by subtraction of the mean
%   res   = res - repmat(mean(res,2),[1 nres 1]);