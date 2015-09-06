function [A,b,x,t] = i_laplace(n,example)
%I_LAPLACE Test problem: inverse Laplace transformation.
%
% [A,b,x,t] = i_laplace(n,example)
%
% Discretization of the inverse Laplace transformation by means of
% Gauss-Laguerre quadrature.  The kernel K is given by
%    K(s,t) = exp(-s*t) ,
% and both integration intervals are [0,inf).
%
% The following examples are implemented, where f denotes
% the solution, and g denotes the right-hand side:
%    1: f(t) = exp(-t/2),        g(s) = 1/(s + 0.5)
%    2: f(t) = 1 - exp(-t/2),    g(s) = 1/s - 1/(s + 0.5)
%    3: f(t) = t^2*exp(-t/2),    g(s) = 2/(s + 0.5)^3
%    4: f(t) = | 0 , t <= 2,     g(s) = exp(-2*s)/s.
%              | 1 , t >  2
%
% The quadrature points are returned in the vector t.

% Reference: J. M. Varah, "Pitfalls in the numerical solution of linear
% ill-posed problems", SIAM J. Sci. Stat. Comput. 4 (1983), 164-176.

% Per Christian Hansen, IMM, Oct. 21, 2006.

% Initialization.
if (n <= 0), error('The order n must be positive'); end
if (nargin == 1), example = 1; end

% Compute equidistand collocation points s.
s = (10/n)*(1:n)';

% Compute abscissas t and weights v from the eigensystem of the
% symmetric tridiagonal system derived from the recurrence
% relation for the Laguerre polynomials.  Sorting of the
% eigenvalues and -vectors is necessary.
t = diag(2*(1:n)-1) - diag((1:n-1),1) - diag((1:n-1),-1);
[Q,t] = eig(t); t = diag(t); [t,indx] = sort(t);
v = abs(Q(1,indx)); clear Q
nz = find(v~=0);

% Set up the coefficient matrix A.  Due to limitations cause
% by finite-precision arithmetic, A has zero columns of n > 195
A = zeros(n,n);
for i=1:n
  for j=nz
    A(i,j) = (1-s(i))*t(j) + 2*log(v(j));
  end
end
A(:,nz) = exp(A(:,nz));

% Compute the right-hand side b and the solution x by means of
% simple collocation.
if (example==1)
  b = ones(n,1)./(s + .5);
  x = exp(-t/2);
elseif (example==2)
  b = ones(n,1)./s - ones(n,1)./(s + .5);
  x = ones(n,1) - exp(-t/2);
elseif (example==3)
  b = 2*ones(n,1)./((s + .5).^3);
  x = (t.^2).*exp(-t/2);
elseif (example==4)
  b = exp(-2*s)./s;
  x = ones(n,1); f = find(t<=2); x(f) = zeros(length(f),1);
else
  error('Illegal example')
end