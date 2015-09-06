function [U,s,V] = csvd(A,tst)
%CSVD Compact singular value decomposition.
%
% s = csvd(A)
% [U,s,V] = csvd(A)
% [U,s,V] = csvd(A,'full')
%
% Computes the compact form of the SVD of A:
%    A = U*diag(s)*V',
% where
%    U  is  m-by-min(m,n)
%    s  is  min(m,n)-by-1
%    V  is  n-by-min(m,n).
%
% If a second argument is present, the full U and V are returned.

% Per Christian Hansen, IMM, 06/22/93.

if (nargin==1)
  if (nargout > 1)
    [m,n] = size(A);
    if (m >= n)
      [U,s,V] = svd(full(A),0); s = diag(s);
    else
      [V,s,U] = svd(full(A)',0); s = diag(s);
    end
  else
    U = svd(full(A));
  end
else
  if (nargout > 1)
    [U,s,V] = svd(full(A)); s = diag(s);
  else
    U = svd(full(A));
  end
end