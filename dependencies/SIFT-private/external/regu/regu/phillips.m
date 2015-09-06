function [A,b,x] = phillips(n)
%PHILLIPS Test problem: Phillips' "famous" problem.
%
% [A,b,x] = phillips(n)
%
% Discretization of the `famous' first-kind Fredholm integral
% equation deviced by D. L. Phillips.  Define the function
%    phi(x) = | 1 + cos(x*pi/3) ,  |x| <  3 .
%             | 0               ,  |x| >= 3
% Then the kernel K, the solution f, and the right-hand side
% g are given by:
%    K(s,t) = phi(s-t) ,
%    f(t)   = phi(t) ,
%    g(s)   = (6-|s|)*(1+.5*cos(s*pi/3)) + 9/(2*pi)*sin(|s|*pi/3) .
% Both integration intervals are [-6,6].
%
% The order n must be a multiple of 4.

% Reference: D. L. Phillips, "A technique for the numerical solution
% of certain integral equations of the first kind", J. ACM 9
% (1962), 84-97.

% Discretized by Galerkin method with orthonormal box functions.

% Per Christian Hansen, IMM, 09/17/92.

% Check input.
if (rem(n,4)~=0), error('The order n must be a multiple of 4'), end

% Compute the matrix A.
h = 12/n; n4 = n/4; r1 = zeros(1,n);
c = cos((-1:n4)*4*pi/n);
r1(1:n4) = h + 9/(h*pi^2)*(2*c(2:n4+1) - c(1:n4) - c(3:n4+2));
r1(n4+1) = h/2 + 9/(h*pi^2)*(cos(4*pi/n)-1);
A = toeplitz(r1);

% Compute the right-hand side b.
if (nargout>1),
  b = zeros(n,1); c = pi/3;
  for i=n/2+1:n
    t1 = -6 + i*h; t2 = t1 - h;
    b(i) =   t1*(6-abs(t1)/2) ...
           + ((3-abs(t1)/2)*sin(c*t1) - 2/c*(cos(c*t1) - 1))/c ...
           - t2*(6-abs(t2)/2) ...
           - ((3-abs(t2)/2)*sin(c*t2) - 2/c*(cos(c*t2) - 1))/c;
    b(n-i+1) = b(i);
  end
  b = b/sqrt(h);
end

% Compute the solution x.
if (nargout==3)
  x = zeros(n,1);
  x(2*n4+1:3*n4) = (h + diff(sin((0:h:(3+10*eps))'*c))/c)/sqrt(h);
  x(n4+1:2*n4) = x(3*n4:-1:2*n4+1);
end