% [z,zu,ldA,Q,T] = DIAGINV_FULL(X,R,B,P)
%
% Computes diag( B*inv(A)*B' ) with A = X'*R*X + B'*P*B.
%
% ARGUMENTS
% X      [m,n]  matrix or operator, m can be zero
% R      [m,m]  square matrix
%     or [m,1]  vector of diagonal          R = diag(R(:))
%     or [1,1]  scalar multiple of identity R = R*eye(n)
%
% B      [q,n]  matrix or operator
% P      [q,q]  square matrix
%     or [q,1]  vector of diagonal          P = diag(P(:))
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
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 21

function [z,zu,ldA,Q,T] = diaginv_full(X,R,B,P)

n = max(size(X,2),size(B,2));

A = zeros(n); ei = zeros(n,1);
if isnumeric(X)
  if numel(R)==size(X,1)
    A = A + X'*(repmat(R(:),1,n).*X);
  else
    A = A + X'*(R*X);
  end
else
  for i=1:n, ei(i)=1; A(:,i) = A(:,i) + mvm(X,R,ei); ei(i)=0; end
end
if isnumeric(B)
  if numel(B)==1
    if numel(P)==size(A,1)
      A = A + (B*B)*diag(P(:));
    elseif numel(P)==1
      A = A + (B*B*P)*eye(n);
    else
      A = A + (B*B)*P;
    end
  elseif numel(P)==size(B,1)
    A = A + B'*(repmat(P(:),1,n).*B);
  else
    A = A + B'*(P*B);
  end
else
  for i=1:n, ei(i)=1; A(:,i) = A(:,i) + mvm(B,P,ei); ei(i)=0; end
end

L = chol((A+A')/2); z = 0; zu = 0;
for i=1:n, ei(i)=1; w=L\ei; zu=zu+w.*conj(w); v=B*w; z=z+v.*conj(v); ei(i)=0;end
if nargout>2, ldA = 2*sum(log(diag(L))); end
if nargout>3, [Q,T] = eig(A); end

function v = mvm(X,R,w)                                        % MVM with X'*R*X
  Xw = X*w; if numel(R)==numel(Xw), RXw = R(:).*Xw; else RXw = R*Xw; end
  v = [X']*RXw;                % the [] around the transpose is needed by Octave