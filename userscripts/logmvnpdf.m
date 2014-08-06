function [v,gM,gS,hMM,hSS,hMS,hSM] = logmvnpdf(X,mu,Sigma,isInv,isChol)
% [Value,GradientM,GradientS,Hessians...] = logmvnpdf(X,mu,Sigma,isInv,isChol)
% Derivatives of the log multivariate normal pdf at points X for mean mu and covariance Sigma.
%
% In:
%   X : points at which to evaluate the log-pdf [#dims x #obs]
%
%   mu : mean vector [#dims x 1]
%
%   S : covariance matrix representation [#dims x #dims]
%
%   isInv : whether Sigma is given in inverted form (e.g., precision matrix) (default: false)
%
%   isChol : whether we have a cholesky factor of Sigma (default: false)
%
% Out:
%   Value : the log-pdf evaluated at X [1 x #obs]
%
%   GradientM : the gradient of the log-pdf with respect to the mean [#dims x #obs]
%
%   GradientS : the gradient of the log-pdf with respect to the covariance matrix 
%               (irrespective of how it's parameterized in the input) [#dims x #dims]
%
%   HessianMM : the Hessian of the log-pdf with respect to the mean [#dims x #dims]
%
%   HessianSS : the Hessian of the log-pdf with respect to the (vectorized) covariance matrix
%               [#dims^2 x #dims^2]
%
%   HessianMS : the partial second derivatives of the log-pdf with respect to 
%               the mean and (vectorized) covariance matrix (off-diagonal block of the Hessian)
%               matrix of the logpdf taken as a function of stacked mean and vectorized covariance
%               [#dims x #dims^2 x #obs]
%
%   HessianSM : transpose of HessianMS
%               [#dims^2 x #dims x #obs]
%

[d,n] = size(X);
if nargin < 2 || isempty(mu)
    mu = zeros(1,d); end
if nargin < 3 || isempty(Sigma)
    Sigma = eye(d,d); end
if nargin < 4 || isempty(isChol)
    isChol = false; end
if nargin < 5 || isempty(isInv)
    isInv = false; end

if ~isvector(mu) && length(mu)==d
    error('The mean should be a vector of length %i but was: %s',d,num2str(size(mu))); end
if ~isequal([d,d],size(Sigma))
    error('The covariance should be of size [%i %i] but was: %s',d,d,num2str(size(S))); end
if isChol && d>1 && any(Sigma(0~=cumsum(eye(d))-eye(d)))
    error('The given cholesky factor should be an upper triangular matrix.'); end

% get Cholesky factor of Sigma
if ~isChol
    C = chol(Sigma); 
else
    C = Sigma; clear Sigma;
end

D = bsxfun(@minus,X,mu(:));
lbase = -0.5*d*log(2*pi);   % log of base measure
ldet = -sum(log(diag(C)));  % log of |Sigma|^-1/2
if isInv
    ldet = -ldet;           % Sigma is negated
    F = C*D;
else
    F = C'\D;
end
v = lbase + ldet -0.5*sum(F.*F);

if nargout>1
    % construct Sigma from its Cholesky factor if necessary
    if isChol
        Sigma = C'*C; end
    % if Sigma is already inverted, the ldivide turns into a product
    if isInv
        gM = Sigma*D;
    else
        gM = Sigma\D;
    end
end
if nargout > 2
    % from here on we need to make sure that we have the actual Sigma and its inverse straight
    if isInv
        iSigma = Sigma;
        Sigma = inv(iSigma);
    end
    gS = -0.5*(Sigma - D*D'); 
end
if nargout > 3
    if ~exist('iSigma','var')
        iSigma = inv(Sigma); end    
    hMM = -iSigma;
end
if nargout > 4
    tmp = -0.5*Sigma;
    hSS = reshape(bsxfun(@times,permute(tmp,[3 1 2 4]),permute(tmp,[2 4 3 1])),d^2,d^2); 
end
if nargout > 5
    A = reshape(bsxfun(@times,eye(d),reshape(D,[1 1 d n])),d,[],n);
    B = reshape(bsxfun(@times,reshape(D,[1 d 1 n]),reshape(eye(d),[d 1 d])),d,[],n);
    hMS = -0.5 * (A+B);
    hSM = permute(hMS,[2 1 3]);
end
