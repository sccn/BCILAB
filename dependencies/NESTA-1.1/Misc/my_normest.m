function [e,cnt] = my_normest(S,St,n,tol, maxiter)
% [norm,iter_count] = my_normest( A, At, n, [tol], [maxiter] )
%   estimates the spectral norm of A using the power method,
%   where A is a function handle to compute A(x) = A*x
%   and At is a function handle to compute At(x) = A'*x
%
% Copied from MATLAB's "normest" function, but allows function handles, not just sparse matrices
    if nargin < 4, tol = 1.e-6; end
    if nargin < 5, maxiter = 20; end
    if isempty(St)
        St = S;  % we assume the matrix is symmetric;
    end
    x = ones(n,1);
    cnt = 0;
    e = norm(x);
    if e == 0, return, end
    x = x/e;
    e0 = 0;
    while abs(e-e0) > tol*e && cnt < maxiter
       e0 = e;
       Sx = S(x);
       if nnz(Sx) == 0
          Sx = rand(size(Sx));
       end
       e = norm(Sx);
       x = St(Sx);
       x = x/norm(x);
       cnt = cnt+1;
    end
end