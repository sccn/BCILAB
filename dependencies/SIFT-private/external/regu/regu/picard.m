function eta = picard(U,s,b,d)
%PICARD Visual inspection of the Picard condition.
%
% eta = picard(U,s,b,d)
% eta = picard(U,sm,b,d)  ,  sm = [sigma,mu]
%
% Plots the singular values, s(i), the abs. value of the Fourier
% coefficients, |U(:,i)'*b|, and a (possibly smoothed) curve of
% the solution coefficients eta(i) = |U(:,i)'*b|/s(i).
%
% If s = [sigma,mu], where gamma = sigma./mu are the generalized
% singular values, then this routine plots gamma(i), |U(:,i)'*b|,
% and (smoothed) eta(i) = |U(:,i)'*b|/gamma(i).
%
% The smoothing is a geometric mean over 2*d+1 points, centered
% at point # i. If nargin = 3, then d = 0 (i.e, no smothing).

% Reference: P. C. Hansen, "The discrete Picard condition for
% discrete ill-posed problems", BIT 30 (1990), 658-672.

% Per Christian Hansen, IMM, April 14, 2001.

% Initialization.
[n,ps] = size(s); beta = abs(U(:,1:n)'*b); eta = zeros(n,1);
if (nargin==3), d = 0; end;
if (ps==2), s = s(:,1)./s(:,2); end
d21 = 2*d+1; keta = 1+d:n-d;
if ~all(s), warning('Division by zero singular values'), end
w = warning('off');
for i=keta
  eta(i) = (prod(beta(i-d:i+d))^(1/d21))/s(i);
end
warning(w);

% Plot the data.
semilogy(1:n,s,'.-',1:n,beta,'x',keta,eta(keta),'o')
xlabel('i')
title('Picard plot')
if (ps==1)
  legend('\sigma_i','|u_i^Tb|','|u_i^Tb|/\sigma_i')
else
  legend('\sigma_i/\mu_i','|u_i^Tb|','|u_i^Tb| / (\sigma_i/\mu_i)')
end