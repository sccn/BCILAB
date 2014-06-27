% c = LINSOLVE_WOODBURY(X,R,P,b)
%
% Computes c = A\b with A = X'*R*X + P positive definite where X is an 
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
% b      [n,1]  vector
%
% RESULTS
% c      [n,1]  vector c = A\b, A = X'*R*X + P
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 11

function c = linsolve_woodbury(X,R,P,b)

[m,n] = size(X);
if m<n                       % use Woodbury identity
  ep = 1e-9*max(P(:));                                % stabilise inversion of P
  if numel(P)==1                       % diag(inv(P))
    c = b/(P+ep);              XiP = X/(P+ep);                   % P is a scalar
  elseif numel(P)==n
    c = b./(P+ep);             XiP = X*sparse(1:n,1:n,1./(P+ep));% P is a vector
  else
    iP = (P+ep*eye(n))\eye(n); XiP = X*iP; c = iP*b;             % P is a matrix
  end
  if numel(P)==1                       % C = inv(R)+X*inv(P)*X'
    C = (X*X')/(P+ep);                                           % P is a scalar
  elseif numel(P)==n
    C = X*diag(1./(P(:)+ep))*X';                                 % P is a vector
  else
    C = X*((P+ep*eye(n))\X');                                    % P is a matrix
  end
  ep = 1e-9*max(R(:));                                % stabilise inversion of R
  if numel(R)==1
    C = C + eye(m)/(R+ep);                                       % R is a scalar
  elseif numel(R)==m
    C = C + diag(1./(R(:)+ep));                                  % R is a vector
  else
    C = C + (R+ep*eye(m))\eye(m);                                % R is a matrix
  end  
  V = XiP'/chol(C);
  c = c - V*(V'*b);
else                         % do direct computation
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
  L = chol(A);
  c = L\(L'\b);
end