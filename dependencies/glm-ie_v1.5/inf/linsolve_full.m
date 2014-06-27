% c = LINSOLVE_FULL(X,R,B,P,b)
%
% Computes c = A\b with A = X'*R*X + B*P*B by
%
% ARGUMENTS
% X      [m,n]  matrix or operator, m can be zero
% R      [m,m]  square matrix
%     or [m,1]  vector of diagonal          R = diag(R(:))
%     or [1,1]  scalar multiple of identity R = R*eye(n)
%
% B      [q,n]  matrix or operator
% P      [n,n]  square matrix
%     or [m,1]  vector of diagonal          P = diag(P(:))
%     or [1,1]  scalar multiple of identity P = P*eye(n)
%
% b      [n,1]  vector
% 
% RESULTS
% c      [n,1]  vector c = A\b, A = X'*X/s2 + B*diag(pi)*B
%
% The matrices X and B can be given implicitely through their multiplication.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 September 25

function c = linsolve_full(X,R,B,P,b)

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
  if numel(P)==size(B,1)
    A = A + B'*(repmat(P(:),1,n).*B);
  else
    A = A + B'*(P*B);
  end
else
  for i=1:n, ei(i)=1; A(:,i) = A(:,i) + mvm(B,P,ei); ei(i)=0; end
end

c = A\b;

function v = mvm(X,R, w)                                       % MVM with X'*R*X
  Xw = X*w; if numel(R)==numel(Xw), RXw = R(:).*Xw; else RXw = R*Xw; end
  v = [X']*RXw;                % the [] around the transpose is needed by Octave