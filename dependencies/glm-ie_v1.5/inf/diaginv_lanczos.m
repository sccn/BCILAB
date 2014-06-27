% [z,zu,ldT,Q,T] = DIAGINV_LANCZOS(X,R,B,P,kmax,output)
%
% Approximates diag( B*inv(A)*B' ) with A = X'*R*X + B'*P*B through
% Krylov subspace estimation via Lanczos iterations using full Gram-Schmidt 
% orthogonalisation in every step.
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
% kmax   [1,1]  number of Lanczos vectors, kmax <= n             [default 100]
% output [1,1]  flag saying whether some output is shown         [default false]
% 
% RESULTS
% z      [q,1]  vector 0 <= z  <= diag(B*inv(A)*B'), A = X'*R*X + B'*P*B
% zu     [n,1]  vector 0 <= zu <= diag(  inv(A)*  ), A = X'*R*X + B'*P*B
% ldT    [1,1]  log determinant of T, ldT <= log(det(A)) if a positive definite
% Q      [n,k]  orthogonal matrix of Lanczos vectors, Q'*Q = I
% T      [k,k]  tridiagonal symmetric matrix, T = Q'*A*Q
%
% The eigensystem of T given by [W,D] = eig(T) can be used to approximate the
% eigensystem of A because [Q*W,D] converges to an eigensystem of A.
% As a result, inv(A) can be approximated by Q*inv(T)*Q'.
%
% The matrices X and B can be given implicitely through their multiplication.
%
% The algorithm is from Schneider & Willsky 2000: Krylov subspace estimation.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 21

function [z,zu,ldT,Q,T] = diaginv_lanczos(X,R,B,P,kmax,output)

n = max(size(X,2),size(B,2));
if nargin<5, kmax = 100; end, kmax = min(n,kmax);              % max. no of MVMs
if nargin<6, output = 0; else zold = 0; end
w = randn(n,1); w = w/norm(w,'fro');                    % random starting vector

dev = @(a,b) norm(a-b)/max([1,norm(a),norm(b)]);            % relative deviation

% $9.2.1 from Golub & van Loan:  only v, w and t are needed as storage of size n
alpha = zeros(1,kmax); beta  = zeros(1,kmax-1);                          % for T
if nargout>1, ldT = 0; end
v = zeros(n,1); b = 1; k = 0;
Q = zeros(n,kmax); z = 0; zu = 0;           % memory for the Lanczos and results
while abs(b)>1e-10 && k<kmax                        % do at most kmax iterations
  if k>0, t = w; w = v/b; v = -b*t; end   
  w(:) = w(:) - Q*(Q'*w(:));                    % Gram-Schmidt orthogonalisation
  v = v + mvmA(X,R,B,P, w);                            % MVM with A: v = v + A*w
  k = k+1;
  a = real(w(:)'*v(:)); alpha(k) = a;                             % needed for T
  v = v - a*w;
  if k==1
    Lkk = sqrt(a);                           % L(k,k), L = chol(T)' <=> L*L' = T
    p   = w/Lkk;                             % Q(:,k) = w
  else
    Lkk1 = b/Lkk;                            % L(k,k-1)
    Lkk  = sqrt( a - Lkk1*Lkk1 );            % L(k,k)
    p    = (w-Lkk1*p)/Lkk;
  end
  if nargout>1 && k>0, ldT = ldT + 2*log(Lkk); end
  b = norm(v,'fro');         if k<kmax, beta(k) = b; end          % needed for T
  Q(:,k) = w(:);                         % store Lanczos vectors for later usage
  zu = zu + p.*conj(p);
  Bp = B*p; z = z + Bp.*conj(Bp);      % filtered backprojected search direction
  if output && ~mod(k,50)
    ndigits = floor(1+log(kmax)/log(10));
    dz = dev(z,zold); zold=z; fprintf('  %*d/%d: dz=%1.4e\n',ndigits,k,kmax,dz)
  end
end
% truncate alpha, beta, Q in case Lanczos did converge earlier
if k<kmax, alpha = alpha(1:k); beta = beta(1:k-1); Q = Q(:,1:k); end
% construct the non-sparse return matrix T
if nargout>4, T = diag(alpha)+diag(beta,+1)+diag(beta,-1); end

function v = mvmA(X,R,B,P, w)             % MVM with A where A = X'*R*X + B'*P*B
  Xw = X*w; if numel(R)==numel(Xw), RXw = R(:).*Xw; else RXw = R*Xw; end
  Bw = B*w; if numel(P)==numel(Bw), PBw = P(:).*Bw; else PBw = P*Bw; end
  v = [X']*RXw + [B']*PBw;     % the [] around the transpose is needed by Octave