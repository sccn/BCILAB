function [U,sm,X,V,W] = cgsvd(A,L)
%CGSVD Compact generalized SVD of a matrix pair in regularization problems.
%
% sm = cgsvd(A,L)
% [U,sm,X,V] = cgsvd(A,L) ,  sm = [sigma,mu]
% [U,sm,X,V,W] = cgsvd(A,L) ,  sm = [sigma,mu]
%
% Computes the generalized SVD of the matrix pair (A,L). The dimensions of
% A and L must be such that [A;L] does not have fewer rows than columns.
%
% If m >= n >= p then the GSVD has the form:
%    [ A ] = [ U  0 ]*[ diag(sigma)      0    ]*inv(X)
%    [ L ]   [ 0  V ] [      0       eye(n-p) ]
%                     [  diag(mu)        0    ]
% where
%    U  is  m-by-n ,    sigma  is  p-by-1
%    V  is  p-by-p ,    mu     is  p-by-1
%    X  is  n-by-n .
%
% Otherwise the GSVD has a more complicated form (see manual for details).
%
% A possible fifth output argument returns W = inv(X).
 
% Reference: C. F. Van Loan, "Computing the CS and the generalized 
% singular value decomposition", Numer. Math. 46 (1985), 479-491. 
 
% Per Christian Hansen, IMM, March 17, 2008. 
 
% Initialization.
[m,n] = size(A); [p,n1] = size(L);
if (n1 ~= n)
  error('No. columns in A and L must be the same')
end
if (m+p < n)
  error('Dimensions must satisfy m+p >= n')
end

% Call Matlab's GSVD routine.
[U,V,W,C,S] = gsvd(full(A),full(L),0);

if (m >= n)
  % The overdetermined or square case.
  sm = [diag(C(1:p,1:p)),diag(S(1:p,1:p))]; 
  if (nargout < 2) 
    U = sm; 
  else 
    % Full decomposition. 
    X = inv(W'); 
  end
else
  % The underdetermined case.
  sm = [diag(C(1:m+p-n,n-m+1:p)),diag(S(n-m+1:p,n-m+1:p))]; 
  if (nargout < 2) 
    U = sm; 
  else 
    % Full decomposition. 
    X = inv(W');
    X = X(:,n-m+1:n); 
  end
end

if (nargout==5), W = W'; end