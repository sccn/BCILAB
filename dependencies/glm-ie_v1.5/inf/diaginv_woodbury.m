% [z,zu,ldA,Q,T] = DIAGINV_WOODBURY(X,R,P)
%
% Computes z = diag(inv(A)) with A = X'*R*X + P positive definite where X is an 
%                     m-by-n matrix and both R and P are positive semi-definite.
%
% We use the Woodbury matrix identity
%   inv(X'*R*X + P) = inv(P) - inv(P)*X'*inv( inv(R)+X*inv(P)*X' )*X*inv(P)
% whenever m<n. Note that we require R and P to be invertible for m<n.
%
% The computational cost is O(N*M^2)>O(M^3) where N=max(m,n) and M=min(m,n)
% due to the Cholesky decomposition on a matrix of size MxM and a matrix-matrix
% multiplication X'*X or X*X'.
%
% ARGUMENTS
% X      [m,n]  full or sparse matrix, no operator allowed
%
% R      [m,m]  square matrix
%     or [m,1]  vector of diagonal          R = diag(R(:))
%     or [1,1]  scalar multiple of identity R = R*eye(n)
%
% P      [n,n]  square matrix
%     or [m,1]  vector of diagonal          P = diag(P(:))
%     or [1,1]  scalar multiple of identity P = P*eye(n)
%
% RESULTS
% z      [q,1]  vector z  = diag(B*inv(A)*B'), A = X'*R*X + B'*P*B, B = I
% zu     [n,1]  vector zu = diag(  inv(A)   ), A = X'*R*X + B'*P*B, B = I
% ldA    [1,1]  log determinant of A
% Q      [n,n]  orthogonal matrix of Lanczos vectors, Q'*Q = I    (for n<m only)
% T      [n,n]  diagonal matrix, T = Q'*A*Q                       (for n<m only)
%
% A is given by Q*T*Q', Q'*Q = I. As a result, inv(A) equals Q*inv(T)*Q'.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 21

function [z,zu,ldA,Q,T] = diaginv_woodbury(X,R,P)

[m,n] = size(X);
if m<n                       % use Woodbury identity
  ep = 1e-9*max(P(:));                               % stabilise inversion of P
  if numel(P)==1                       % diag(inv(P))
    z = 1/(P+ep);              XiP = X*z;                        % P is a scalar
  elseif numel(P)==n
    z = 1./(P+ep);             XiP = X*sparse(1:n,1:n,z);        % P is a vector
  else
    iP = (P+ep*eye(n))\eye(n); XiP = X*iP; z = diag(iP);         % P is a matrix
  end
  if numel(P)==1                       % C = inv(R)+X*inv(P)*X'
    C = (X*X')/(P+ep);                                           % P is a scalar
    if nargout>1, ldP = n*log(P+ep); end
  elseif numel(P)==n
    C = X*diag(1./(P(:)+ep))*X';                                 % P is a vector
    if nargout>1, ldP = sum(log(P(:)+ep)); end
  else
    C = X*((P+ep*eye(n))\X');                                    % P is a matrix
    if nargout>1, L = chol(P+ep*eye(n)); ldP = 2*sum(log(diag(L))); end
  end
  ep = 1e-9*max(R(:));                                % stabilise inversion of R
  if numel(R)==1
    C = C + eye(m)/(R+ep);                                       % R is a scalar
    if nargout>1, ldR = m*log(R+ep); end
  elseif numel(R)==m
    C = C + diag(1./(R(:)+ep));                                  % R is a vector
    if nargout>1, ldR = sum(log(R(:)+ep)); end
  else
    C = C + (R+ep*eye(m))\eye(m);                                % R is a matrix
    if nargout>1, L = chol(R+ep*eye(m)); ldR = 2*sum(log(diag(L))); end
  end
  L = chol(C); V = XiP'/L; z = z - sum(V.*V,2);
  if nargout>1, zu = z; end
  if nargout>2, ldA = ldP + ldR + 2*sum(log(diag(L))); end
  Q = []; T = [];                                 % there is no canonical output
else                                                     % do direct computation
  if numel(R)==m
    A = X'*diag(R(:))*X;                                         % R is a vector
  else
    A = X'*R*X;                                      % R is a scalar or a matrix
  end
  if numel(P)==1
    A = A + P*eye(n);                                            % P is a scalar
  elseif numel(P)==n
    A = A + diag(P(:));                                          % P is a vector
  else
    A = A + P;                                                   % P is a matrix
  end
  L = chol(A); iL = L \ eye(n);
  z = sum(iL.*iL,2);
  if nargout>1, zu = z; end
  if nargout>2, ldA = 2*sum(log(diag(L))); end
  if nargout>3, [Q,T] = eig(A); end
end