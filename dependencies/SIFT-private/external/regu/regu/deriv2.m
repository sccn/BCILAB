function [A,b,x] = deriv2(n,example)
%DERIV2 Test problem: computation of the second derivative.
%
% [A,b,x] = deriv2(n,example)
%
% This is a mildly ill-posed problem.  It is a discretization of a
% first kind Fredholm integral equation whose kernel K is the
% Green's function for the second derivative:
%    K(s,t) = | s(t-1)  ,  s <  t .
%             | t(s-1)  ,  s >= t
% Both integration intervals are [0,1], and as right-hand side g
% and correspond solution f one can choose between the following:
%    example = 1 : g(s) = (s^3 - s)/6          ,  f(t) = t
%    example = 2 : g(s) = exp(s) + (1-e)s - 1  ,  f(t) = exp(t)
%    example = 3 : g(s) = | (4s^3 - 3s)/24               ,  s <  0.5
%                         | (-4s^3 + 12s^2 - 9s + 1)/24  ,  s >= 0.5
%                  f(t) = | t    ,  t <  0.5
%                         | 1-t  ,  t >= 0.5

% References.  The first two examples are from L. M. Delves & J. L.
% Mohamed, "Computational Methods for Integral Equations", Cambridge
% University Press, 1985; p. 310.  The third example is from A. K.
% Louis and P. Maass, "A mollifier method for linear operator equations
% of the first kind", Inverse Problems 6 (1990), 427-440.

% Discretized by the Galerkin method with orthonormal box functions.

% Per Christian Hansen, IMM, 04/21/97.

% Initialization.
if (nargin==1), example = 1; end
h = 1/n; sqh = sqrt(h); h32 = h*sqh; h2 = h^2; sqhi = 1/sqh;
t = 2/3; A = zeros(n,n);

% Compute the matrix A.
for i=1:n
  A(i,i) = h2*((i^2 - i + 0.25)*h - (i - t));
  for j=1:i-1
    A(i,j) = h2*(j-0.5)*((i-0.5)*h-1);
  end
end
A = A + tril(A,-1)';

% Compute the right-hand side vector b.
if (nargout>1)
  b = zeros(n,1);
  if (example==1)
    for i=1:n
      b(i) = h32*(i-0.5)*((i^2 + (i-1)^2)*h2/2 - 1)/6;
    end
  elseif (example==2)
    ee = 1 - exp(1);
    for i=1:n
      b(i) = sqhi*(exp(i*h) - exp((i-1)*h) + ee*(i-0.5)*h2 - h);
    end
  elseif (example==3)
    if (rem(n,2)~=0), error('Order n must be even'), else
      for i=1:n/2
        s12 = (i*h)^2; s22 = ((i-1)*h)^2;
        b(i) = sqhi*(s12 + s22 - 1.5)*(s12 - s22)/24;
      end
      for i=n/2+1:n
        s1 = i*h; s12 = s1^2; s2 = (i-1)*h; s22 = s2^2;
        b(i) = sqhi*(-(s12+s22)*(s12-s22) + 4*(s1^3 - s2^3) - ...
                    4.5*(s12 - s22) + h)/24;
      end
    end
  else
    error('Illegal value of example')
  end
end

% Compute the solution vector x.
if (nargout==3)
  x = zeros(n,1);
  if (example==1)
    for i=1:n, x(i) = h32*(i-0.5); end
  elseif(example==2)
    for i=1:n, x(i) = sqhi*(exp(i*h) - exp((i-1)*h)); end
  else
    for i=1:n/2,   x(i) = sqhi*((i*h)^2 - ((i-1)*h)^2)/2; end
    for i=n/2+1:n, x(i) = sqhi*(h - ((i*h)^2 - ((i-1)*h)^2)/2); end
  end
end