% c = LINSOLVE_LCG(X,R,B,P,b,kmax,c0)
%
% Computes c = A\b with A = X'*R*X + B'*P*B by Linear Conjugate Gradients (LCG) 
% as described in § 10.2.6 of Golub, van Loan, Matrix Computations with at most
% kmax MVMs starting from an intial value c0.
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
% kmax   [1,1]  scalar, default value is n
% c0     [n,1]  vector, default is zeros(n,1)
% 
% RESULTS
% c      [n,1]  vector c = A\b, A = X'*R*X + B'*P*B
% 
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 13

function [c,k,res] = linsolve_lcg(X,R,B,P,b,kmax,c0)

if nargin>6, c=c0; else c=zeros(size(b)); end
if nargin<6, kmax=length(b); end 

r  = b - mvmA(X,R,B,P, c);                 % r = b-A*c;
nb = norm(b,'fro');                        % |b|
k  = 0;                                    % iteration counter

while norm(r,'fro')>eps*nb && k<kmax
  k = k+1;
  if k==1
    p = -r;
  else
    beta = (r(:)'*w(:))/(p(:)'*w(:));     
    p = -r + beta*p;
  end
  w = mvmA(X,R,B,P, p);    % w = A*p;
  alpha = (r(:)'*p(:))/(p(:)'*w(:));
  c = c + alpha*p;
  r = r - alpha*w;
end
res = norm(r,'fro')/nb;                    % residual norm divided by norm of b

function v = mvmA(X,R,B,P, w)             % MVM with A where A = X'*R*X + B'*P*B
  Xw = X*w; if numel(R)==numel(Xw), RXw = R(:).*Xw; else RXw = R*Xw; end
  Bw = B*w; if numel(P)==numel(Bw), PBw = P(:).*Bw; else PBw = P*Bw; end
  v = [X']*RXw + [B']*PBw;     % the [] around the transpose is needed by Octave