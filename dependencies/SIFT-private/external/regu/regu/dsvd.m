function [x_lambda,rho,eta] = dsvd(U,s,V,b,lambda)
%DSVD Damped SVD and GSVD regularization.
%
% [x_lambda,rho,eta] = dsvd(U,s,V,b,lambda)
% [x_lambda,rho,eta] = dsvd(U,sm,X,b,lambda) ,  sm = [sigma,mu]
%
% Computes the damped SVD solution defined as
%    x_lambda = V*inv(diag(s + lambda))*U'*b .
% If lambda is a vector, then x_lambda is a matrix such that
%    x_lambda = [ x_lambda(1), x_lambda(2), ... ] .
%
% If sm and X are specified, then the damped GSVD solution:
%    x_lambda = X*[ inv(diag(sigma + lambda*mu)) 0 ]*U'*b
%                 [            0                 I ]
% is computed.
%
% The solution norm (standard-form case) or seminorm (general-form
% case) and the residual norm are returned in eta and rho.

% Reference: M. P. Ekstrom & R. L. Rhoads, "On the application of
% eigenvector expansions to numerical deconvolution", J. Comp.
% Phys. 14 (1974), 319-340.
% The extension to GSVD is by P. C. Hansen.

% Per Christian Hansen, IMM, April 14, 2003.

% Initialization.
if (min(lambda)<0)
  error('Illegal regularization parameter lambda')
end
m = size(U,1);
n = size(V,1);
[p,ps] = size(s);
beta = U(:,1:p)'*b;
ll = length(lambda); x_lambda = zeros(n,ll);
rho = zeros(ll,1); eta = zeros(ll,1);

% Treat each lambda separately.
if (ps==1)

  % The standard-form case.
  for i=1:ll
    x_lambda(:,i) = V(:,1:p)*(beta./(s + lambda(i)));
    rho(i) = lambda(i)*norm(beta./(s + lambda(i)));
    eta(i) = norm(x_lambda(:,i));
  end
  if (nargout > 1 && size(U,1) > p)
    rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
  end

elseif (m>=n)

  % The overdetermined or square general-form case.
  x0 = V(:,p+1:n)*U(:,p+1:n)'*b;
  for i=1:ll
    xi = beta./(s(:,1) + lambda(i)*s(:,2));
    x_lambda(:,i) = V(:,1:p)*xi + x0;
    rho(i) = lambda(i)*norm(beta./(s(:,1)./s(:,2) + lambda(i)));
    eta(i) = norm(s(:,2).*xi);
  end
  if (nargout > 1 && size(U,1) > p)
    rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
  end

else

  % The underdetermined general-form case.
  x0 = V(:,p+1:m)*U(:,p+1:m)'*b;
  for i=1:ll
    xi = beta./(s(:,1) + lambda(i)*s(:,2));
    x_lambda(:,i) = V(:,1:p)*xi + x0;
    rho(i) = lambda(i)*norm(beta./(s(:,1)./s(:,2) + lambda(i)));
    eta(i) = norm(s(:,2).*xi);
  end

end