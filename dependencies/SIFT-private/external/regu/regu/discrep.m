function [x_delta,lambda] = discrep(U,s,V,b,delta,x_0)
%DISCREP Discrepancy principle criterion for choosing the reg. parameter.
%
% [x_delta,lambda] = discrep(U,s,V,b,delta,x_0)
% [x_delta,lambda] = discrep(U,sm,X,b,delta,x_0)  ,  sm = [sigma,mu]
%
% Least squares minimization with a quadratic inequality constraint:
%    min || x - x_0 ||       subject to   || A x - b || <= delta
%    min || L (x - x_0) ||   subject to   || A x - b || <= delta
% where x_0 is an initial guess of the solution, and delta is a
% positive constant.  Requires either the compact SVD of A saved as
% U, s, and V, or part of the GSVD of (A,L) saved as U, sm, and X.
% The regularization parameter lambda is also returned.
%
% If delta is a vector, then x_delta is a matrix such that
%    x_delta = [ x_delta(1), x_delta(2), ... ] .
%
% If x_0 is not specified, x_0 = 0 is used.

% Reference: V. A. Morozov, "Methods for Solving Incorrectly Posed
% Problems", Springer, 1984; Chapter 26.

% Per Christian Hansen, IMM, August 6, 2007.

% Initialization.
m = size(U,1);          n = size(V,1);
[p,ps] = size(s);       ld  = length(delta);
x_delta = zeros(n,ld);  lambda = zeros(ld,1);  rho = zeros(p,1);
if (min(delta)<0)
  error('Illegal inequality constraint delta')
end
if (nargin==5), x_0 = zeros(n,1); end
if (ps == 1), omega = V'*x_0; else omega = V\x_0; end

% Compute residual norms corresponding to TSVD/TGSVD.
beta = U'*b;
if (ps == 1)
  delta_0 = norm(b - U*beta);
  rho(p) = delta_0^2;
  for i=p:-1:2
    rho(i-1) = rho(i) + (beta(i) - s(i)*omega(i))^2;
  end
else
  delta_0 = norm(b - U*beta);
  rho(1) = delta_0^2;
  for i=1:p-1
    rho(i+1) = rho(i) + (beta(i) - s(i,1)*omega(i))^2;
  end
end

% Check input.
if (min(delta) < delta_0)
  error('Irrelevant delta < || (I - U*U'')*b ||')
end

% Determine the initial guess via rho-vector, then solve the nonlinear
% equation || b - A x ||^2 - delta_0^2 = 0 via Newton's method.
if (ps == 1)
    
  % The standard-form case.
  s2 = s.^2;
  for k=1:ld
    if (delta(k)^2 >= norm(beta - s.*omega)^2 + delta_0^2)
      x_delta(:,k) = x_0;
    else
      [dummy,kmin] = min(abs(rho - delta(k)^2));
      lambda_0 = s(kmin);
      lambda(k) = newton(lambda_0,delta(k),s,beta,omega,delta_0);
      e = s./(s2 + lambda(k)^2); f = s.*e;
      x_delta(:,k) = V(:,1:p)*(e.*beta + (1-f).*omega);
    end
  end
  
elseif (m>=n)
    
  % The overdetermined or square genera-form case.
  omega = omega(1:p); gamma = s(:,1)./s(:,2);
  x_u   = V(:,p+1:n)*beta(p+1:n);
  for k=1:ld
    if (delta(k)^2 >= norm(beta(1:p) - s(:,1).*omega)^2 + delta_0^2)
      x_delta(:,k) = V*[omega;U(:,p+1:n)'*b];
    else
      [dummy,kmin] = min(abs(rho - delta(k)^2));
      lambda_0 = gamma(kmin);
      lambda(k) = newton(lambda_0,delta(k),s,beta(1:p),omega,delta_0);
      e = gamma./(gamma.^2 + lambda(k)^2); f = gamma.*e;
      x_delta(:,k) = V(:,1:p)*(e.*beta(1:p)./s(:,2) + ...
                               (1-f).*s(:,2).*omega) + x_u;
    end
  end
  
else

  % The underdetermined general-form case.
  omega = omega(1:p); gamma = s(:,1)./s(:,2);
  x_u   = V(:,p+1:m)*beta(p+1:m);
  for k=1:ld
    if (delta(k)^2 >= norm(beta(1:p) - s(:,1).*omega)^2 + delta_0^2)
      x_delta(:,k) = V*[omega;U(:,p+1:m)'*b];
    else
      [dummy,kmin] = min(abs(rho - delta(k)^2));
      lambda_0 = gamma(kmin);
      lambda(k) = newton(lambda_0,delta(k),s,beta(1:p),omega,delta_0);
      e = gamma./(gamma.^2 + lambda(k)^2); f = gamma.*e;
      x_delta(:,k) = V(:,1:p)*(e.*beta(1:p)./s(:,2) + ...
                               (1-f).*s(:,2).*omega) + x_u;
    end
  end
end

%-------------------------------------------------------------------

function lambda = newton(lambda_0,delta,s,beta,omega,delta_0)
%NEWTON Newton iteration (utility routine for DISCREP).
%
% lambda = newton(lambda_0,delta,s,beta,omega,delta_0)
%
% Uses Newton iteration to find the solution lambda to the equation
%    || A x_lambda - b || = delta ,
% where x_lambda is the solution defined by Tikhonov regularization.
%
% The initial guess is lambda_0.
%
% The norm || A x_lambda - b || is computed via s, beta, omega and
% delta_0.  Here, s holds either the singular values of A, if L = I,
% or the c,s-pairs of the GSVD of (A,L), if L ~= I.  Moreover,
% beta = U'*b and omega is either V'*x_0 or the first p elements of
% inv(X)*x_0.  Finally, delta_0 is the incompatibility measure.

% Reference: V. A. Morozov, "Methods for Solving Incorrectly Posed
% Problems", Springer, 1984; Chapter 26.

% Per Christian Hansen, IMM, 12/29/97.

% Set defaults.
thr = sqrt(eps);  % Relative stopping criterion.
it_max = 50;      % Max number of iterations.

% Initialization.
if (lambda_0 < 0)
  error('Initial guess lambda_0 must be nonnegative')
end
[p,ps] = size(s);
if (ps==2), sigma = s(:,1); s = s(:,1)./s(:,2); end
s2 = s.^2;

% Use Newton's method to solve || b - A x ||^2 - delta^2 = 0.
% It was found experimentally, that this formulation is superior
% to the formulation || b - A x ||^(-2) - delta^(-2) = 0.
lambda = lambda_0; step = 1; it = 0;
while (abs(step) > thr*lambda & abs(step) > thr & it < it_max), it = it+1;
  f = s2./(s2 + lambda^2);
  if (ps==1)
    r = (1-f).*(beta - s.*omega);
    z = f.*r;
  else
    r = (1-f).*(beta - sigma.*omega);
    z = f.*r;
  end
  step = (lambda/4)*(r'*r + (delta_0+delta)*(delta_0-delta))/(z'*r);
  lambda = lambda - step;
  % If lambda < 0 then restart with smaller initial guess.
  if (lambda < 0), lambda = 0.5*lambda_0; lambda_0 = 0.5*lambda_0; end
end

% Terminate with an error if too many iterations.
if (abs(step) > thr*lambda & abs(step) > thr)
  error(['Max. number of iterations (',num2str(it_max),') reached'])
end