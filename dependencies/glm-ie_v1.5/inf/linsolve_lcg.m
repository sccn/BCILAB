% [c,k,res,mvmMfun,prec] = LINSOLVE_LCG(X,R,B,P,b,kmax,c,tol,prec)
%
% Computes c = A\b with A = X'*R*X + B'*P*B by Linear Conjugate Gradients (LCG) 
% as described in § 10.2.6 of Golub, van Loan, Matrix Computations with at most
% kmax MVMs.
%
% ARGUMENTS
% X      [m,n]   matrix or operator, m can be zero
% R      [m,m]   square matrix
%     or [m,1]   vector of diagonal          R = diag(R(:))
%     or [1,1]   scalar multiple of identity R = R*eye(n)
%
% B      [q,n]   matrix or operator
% P      [n,n]   square matrix
%     or [m,1]   vector of diagonal          P = diag(P(:))
%     or [1,1]   scalar multiple of identity P = P*eye(n)
%
% b      [n,1]   vector
% kmax   [1,1]   scalar, default value is n
% c      [n,1]   vector, default is zeros(n,1)
% tol    [1,1]   scalar, residual relative tolerance |r|/|b| < tol
% prec   [1,1]   flag, try and use preconditioner
% 
% RESULTS
% c      [n,1]   vector c = A\b, A = X'*R*X + B'*P*B
% k      [1,1]   number of commited iterations
% res    [k+1,1] residual norm
% prec   [1,1]   flag whether preconditioning is possible at all
% 
% (c) by Hannes Nickisch, MPI for Biological Cybernetics and 
%        George Papandreou, University of California, LA, 2011 October 20

function [c,k,res,mvmMfun,prec] = linsolve_lcg(X,R,B,P,b,kmax,c,tol,prec)

if nargin<6, kmax=length(b); end 
if nargin<7 || isempty(c), c=zeros(size(b)); end
if nargin<8 || isempty(tol), tol=1e-10; end
if nargin<9, prec=[]; end

% If |b|=0, no need to iterate
nb = norm(b,'fro');                        % |b|
if nb==0
  c = zeros(size(b));
  k = 0;
  res = 0;
  return
end

% Check if fast Jacobi-DFT preconditioner is possible.
% If it is, then we use as preconditioner the circulant matrix
% M = mean(R)*X'*X + mean(P)*B'*B
% Then M\r is implemented as pointwise division in the DFT domain.
if isempty(prec) || prec
  fftdiag = ~isnumeric(X) && ~isnumeric(B);
  if fftdiag
    dX = diagFAtAFt(X,'cheap'); dB = diagFAtAFt(B,'cheap');
    fftdiag = numel(dX)>0 && numel(dB)>0;
  end
  if fftdiag
    R_mean = meandiag(R);
    P_mean = meandiag(P);
    d = R_mean*dX + P_mean*dB;
    F = matFFTNmask(true(size(d))); d = d(:); % DFT matrix
    prec = 1;
  else
    if prec
      error('linsolve_lcg: requested preconditioning is not possible');
    end
    prec = 0;
  end
end

res = zeros(kmax+1,1);

mvmAfun = @(c) mvmA(X,R,B,P, c);
if prec
  mvmMinvfun = @(r) mvmMinv(F,d, r);
end

% See Shewchuk, pp. 50 and 51 for the algorithm
r = b - mvmAfun(c);                       % r = b-A*c;
res(1) = norm(r,'fro')^2;

if ~prec                                        % proceed without preconditioner
  z  = r;
  for k=1:kmax
    w = mvmAfun(z);
    alpha = res(k)/(z(:)'*w(:));
    c = c + alpha*z;
    r = r - alpha*w;
    res(k+1) = norm(r,'fro')^2;
    beta = res(k+1)/res(k);
    z = r + beta*z;
    if sqrt(res(k+1))<tol*nb
      break
    end
  end
else                                                        % use preconditioner
  z = mvmMinvfun(r);
  rs = r(:)'*z(:);

  for k=1:kmax
    w = mvmAfun(z);
    alpha = rs/(z(:)'*w(:));
    c = c + alpha*z;
    r = r - alpha*w;
    res(k+1) = norm(r,'fro')^2;
    s = mvmMinvfun(r);
    rs_old = rs;
    rs = r(:)'*s(:);
    beta = rs/rs_old;
    z = s + beta*z;
    if sqrt(res(k+1))<tol*nb
      break
    end
  end
  
end

res = sqrt(res(1:k+1));

% also return the preconditioning operator if requested
if nargout>=4
  if ~prec
    mvmMfun = @(x) mvmI(x);
  else
    mvmMfun = @(x) mvmM(F,d,x);
  end
end


function v = mvmA(X,R,B,P, w)             % MVM with A where A = X'*R*X + B'*P*B
  Xw = X*w; if numel(R)==numel(Xw), RXw = R(:).*Xw; else RXw = R*Xw; end
  Bw = B*w; if numel(P)==numel(Bw), PBw = P(:).*Bw; else PBw = P*Bw; end
  v = [X']*RXw + [B']*PBw;     % the [] around the transpose is needed by Octave
  
function z = mvmMinv(F,d, r)                             % M\r in the DFT domain
  if numel(r)>numel(d)
    z = cx2re([F']*((F*re2cx(r))./d));
  else
    z = [F']*((F*r)./d);
  end

function [z,d] = mvmM(F,d, x)                            % M*x in the DFT domain
  if numel(x)>numel(d)
    z = cx2re([F']*((F*re2cx(x)).*d));
  else
    z = [F']*((F*x).*d);
  end

function [x,d] = mvmI(x)                                          % M*x when M=I
  d = 1;  

function J_mean = meandiag(J)                            % Compute mean(diag(J))
  [m n] = size(J);
  if isempty(J)
    J_mean = [];
  elseif m==1 || n==1,
    J_mean = mean(J);
  elseif m==n
    J_mean = mean(diag(J));
  else
    error('meandiag: wrong J');
  end
