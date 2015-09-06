function [X,rho,eta] = art(A,b,k)
%ART  Algebraic reconstruction technique (Kaczmarz's method).
%
% [X,rho,eta] = art(A,b,k)
%
% Classical Kaczmarz iteration, or ART (algebraic reconstruction
% technique), applied to the system A x = b.  The number of
% iterations is k.

% Reference: F. Natterer and F. Wübbeling, Mathematical Methods
% in Image Reconstruction, SIAM, Philadelphia, 2001; Sect. 5.3.1.

% Per Christian Hansen, IMM, Dec. 6, 2006.

% Initialization.
if (k < 1), error('Number of steps k must be positive'), end
[m,n] = size(A); X = zeros(n,k);
if (nargout > 1)
   eta = zeros(k,1); rho = eta;
end

% Prepare for iteration.
x = zeros(n,1);
%nai2 = full(sum(A.*A,2));
nai2 = full(sum(abs(A.*A),2));
I = find(nai2>0)';

% Iterate.
for j=1:k
   for i=I
      Ai = full(A(i,:));
      x = x + (b(i)-Ai*x)*Ai'/nai2(i);
   end
   if (nargout > 1)
      eta(j) = norm(x); rho(j) = norm(b-A*x);
   end
   X(:,j) = x;
end