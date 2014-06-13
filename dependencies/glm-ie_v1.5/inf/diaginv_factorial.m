% [z,zu,ldA,Q,T] = DIAGINV_FULL(X,R,B,P)
%
% Approximates diag( B*inv(A)*B' ) with A = X'*R*X + B'*P*B by making a diagonal
% approximation to A: 
%
% ARGUMENTS
% X      [m,n]  matrix or operator, m can be zero
% R      [m,1]  vector of diagonal          R = diag(R(:))
%     or [1,1]  scalar multiple of identity R = R*eye(n)
%
% B      [q,n]  matrix or operator
% P      [q,1]  vector of diagonal          P = diag(P(:))
%     or [1,1]  scalar multiple of identity P = P*eye(n)
% 
% RESULTS
% z      [q,1]  vector z  = diag(B*inv(A)*B'), A = X'*R*X + B'*P*B
% zu     [n,1]  vector zu = diag(  inv(A)   ), A = X'*R*X + B'*P*B
% ldA    [1,1]  log determinant of A
% Q      [n,n]  orthogonal matrix of Lanczos vectors, Q'*Q = I
% T      [n,n]  diagonal matrix, T = Q'*A*Q
%
% A is given by Q*T*Q', Q'*Q = I. As a result, inv(A) equals Q*inv(T)*Q'.
%
% The matrices X and B can be given implicitely through their multiplication.
%
% (c) by Hannes Nickisch, Philips Research, 2012 January 13

function [z,zu,ldA,Q,T] = diaginv_factorial(X,R,B,P)

n = max(size(X,2),size(B,2));
m = size(X,1); q = size(B,1);

dgA = zeros(n,1);
if isnumeric(X)
  if numel(R)==m
    dgA = dgA + diag( X'*(repmat(R(:),1,n).*X) );
  else
    dgA = dgA + diag( X'*(R*X) );
  end
else
  dgA = dgA + reshape(mvmAA(X,R.*ones(m,1),true),[],1);
end
if isnumeric(B)
  if numel(B)==1
    if numel(P)==n
      dgA = dgA +  diag( (B*B)*diag(P(:)) );
    elseif numel(P)==1
      dgA = dgA +  diag( (B*B*P)*eye(n) );
    else
      dgA = dgA +  diag( (B*B)*P );
    end
  elseif numel(P)==q
    dgA = dgA +  diag( B'*(repmat(P(:),1,n).*B) );
  else
    dgA = dgA +  diag( B'*(P*B) );
  end
else
  dgA = dgA + reshape(mvmAA(B,P.*ones(q,1),true),[],1);
end
zu = 1./dgA;
if nargout>1
  if isnumeric(B)
    z = diag(B*diag(zu)*B');
  else
    z = mvmAA(B,zu);
  end
end
if nargout>2, ldA = sum(log(dgA)); end
if nargout>3, Q = speye(n); end
if nargout>4, T = sparse(1:n,1:n,dgA); end
