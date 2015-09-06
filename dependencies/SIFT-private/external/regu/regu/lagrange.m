function [La,dLa,lambda0] = lagrange(U,s,b,more)
%LAGRANGE Plot the Lagrange function for Tikhonov regularization.
%
% [La,dLa,lambda0] = lagrange(U,s,b,more)
% [La,dLa,lambda0] = lagrange(U,sm,b,more)  ,  sm = [sigma,mu]
%
% Plots the Lagrange function
%    La(lambda) = || A x - b ||^2 + lambda^2*|| L x ||^2
% and its first derivative dLa = dLa/dlambda versus lambda.
% Here, x is the Tikhonov regularized solution.
%
% If nargin = 4, || A x - b || and || L x || are also plotted.
%
% Returns La, dLa, and the value lambda0 of lambda for which
% dLa has its minimum.

% Per Christian Hansen, IMM, Feb. 21, 2001.

% Set default number of points.
npoints = 200;

% Initialization.
[m,n] = size(U); [p,ps] = size(s);
beta = U'*b; beta2 = norm(b)^2 - norm(beta)^2;
if (ps==2)
  s = s(p:-1:1,1)./s(p:-1:1,2); beta = beta(p:-1:1);
end
xi = beta(1:p)./s;

% Compute the L-curve.
eta = zeros(npoints,1); rho = eta;
lambda(npoints,1) = s(p);
ratio = (s(1)/s(p))^(1/(npoints-1));
for i=npoints-1:-1:1, lambda(i) = ratio*lambda(i+1); end
for i=1:npoints
  f = fil_fac(s,lambda(i));
  eta(i) = norm(f.*xi);
  rho(i) = norm((1-f).*beta(1:p));
end
if (m > n && beta2 > 0), rho = sqrt(rho.^2 + beta2); end

% Compute the Lagrange function and its derivative.
La = rho.^2 + (lambda.^2).*(eta.^2);
dLa = 2*lambda.*(eta.^2);
[mindLa,mindLi] = min(dLa); lambda0 = lambda(mindLi);

% Plot the functions.
if (nargin==3)
  loglog(lambda,La,'-',lambda,dLa,'--',lambda0,mindLa,'o')
  legend('La','dLa/d\lambda')
else
  loglog(lambda,La,'-',lambda,dLa,'--',lambda,eta,':',lambda,rho,'-.',...
         lambda0,mindLa,'o')
  legend('La','dLa/d\lambda','|| L x ||_2','|| A x - b ||_2')
end
xlabel('\lambda')
title('Lagrange function La and its derivative')