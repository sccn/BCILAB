function [x_lambda,rho,eta] = tikhonov(U,s,V,b,lambda,x_0)
%TIKHONOV Tikhonov regularization.
%
% [x_lambda,rho,eta] = tikhonov(U,s,V,b,lambda,x_0)
% [x_lambda,rho,eta] = tikhonov(U,sm,X,b,lambda,x_0) ,  sm = [sigma,mu]
%
% Computes the Tikhonov regularized solution x_lambda.  If the SVD
% is used, i.e. if U, s, and V are specified, then standard-form
% regularization is applied:
%    min { || A x - b ||^2 + lambda^2 || x - x_0 ||^2 } .
% If, on the other hand, the GSVD is used, i.e. if U, sm, and X are
% specified, then general-form regularization is applied:
%    min { || A x - b ||^2 + lambda^2 || L (x - x_0) ||^2 } .
%
% If x_0 is not specified, then x_0 = 0 is used
%
% Note that x_0 cannot be used if A is underdetermined and L ~= I.
%
% If lambda is a vector, then x_lambda is a matrix such that
%    x_lambda = [ x_lambda(1), x_lambda(2), ... ] .
%
% The solution norm (standard-form case) or seminorm (general-form
% case) and the residual norm are returned in eta and rho.

% Per Christian Hansen, IMM, April 14, 2003.

% Reference: A. N. Tikhonov & V. Y. Arsenin, "Solutions of
% Ill-Posed Problems", Wiley, 1977.

% Initialization.
if (min(lambda)<0)
  error('Illegal regularization parameter lambda')
end
m = size(U,1);
n = size(V,1);
[p,ps] = size(s);
beta = U(:,1:p)'*b;
zeta = s(:,1).*beta;
ll = length(lambda); x_lambda = zeros(n,ll);
rho = zeros(ll,1); eta = zeros(ll,1);

% Treat each lambda separately.
if (ps==1)

  % The standard-form case.
  if (nargin==6), omega = V'*x_0; end
  for i=1:ll
    if (nargin==5)
      x_lambda(:,i) = V(:,1:p)*(zeta./(s.^2 + lambda(i)^2));
      rho(i) = lambda(i)^2*norm(beta./(s.^2 + lambda(i)^2));
    else
      x_lambda(:,i) = V(:,1:p)*...
        ((zeta + lambda(i)^2*omega)./(s.^2 + lambda(i)^2));
      rho(i) = lambda(i)^2*norm((beta - s.*omega)./(s.^2 + lambda(i)^2));
    end
    eta(i) = norm(x_lambda(:,i));
  end
  if (nargout > 1 && size(U,1) > p)
    rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
  end

elseif (m>=n)

  % The overdetermined or square general-form case.
  gamma2 = (s(:,1)./s(:,2)).^2;
  if (nargin==6), omega = V\x_0; omega = omega(1:p); end
  if (p==n)
    x0 = zeros(n,1);
  else
    x0 = V(:,p+1:n)*U(:,p+1:n)'*b;
  end
  for i=1:ll
    if (nargin==5)
      xi = zeta./(s(:,1).^2 + lambda(i)^2*s(:,2).^2);
      x_lambda(:,i) = V(:,1:p)*xi + x0;
      rho(i) = lambda(i)^2*norm(beta./(gamma2 + lambda(i)^2));
    else
      xi = (zeta + lambda(i)^2*(s(:,2).^2).*omega)./...
           (s(:,1).^2 + lambda(i)^2*s(:,2).^2);
      x_lambda(:,i) = V(:,1:p)*xi + x0;
      rho(i) = lambda(i)^2*norm((beta - s(:,1).*omega)./...
               (gamma2 + lambda(i)^2));
    end
    eta(i) = norm(s(:,2).*xi);
  end
  if (nargout > 1 && size(U,1) > p)
    rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
  end

else

  % The underdetermined general-form case.
  gamma2 = (s(:,1)./s(:,2)).^2;
  if (nargin==6), error('x_0 not allowed'), end
  if (p==m)
    x0 = zeros(n,1);
  else
    x0 = V(:,p+1:m)*U(:,p+1:m)'*b;
  end
  for i=1:ll
    xi = zeta./(s(:,1).^2 + lambda(i)^2*s(:,2).^2);
    x_lambda(:,i) = V(:,1:p)*xi + x0;
    rho(i) = lambda(i)^2*norm(beta./(gamma2 + lambda(i)^2));
    eta(i) = norm(s(:,2).*xi);
  end

end