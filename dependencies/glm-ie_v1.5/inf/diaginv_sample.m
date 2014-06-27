% [z,zu,ldA,Q,T] = DIAGINV_SAMPLE(X,R,B,P,Ns,Ncg)
%
% Approximates diag( B*inv(A)*B' ) with A = X'*R*X + B'*P*B as
% Monte-Carlo sample-average. Note that R and P correspond to inverse
% covariance matrices.
%
% Specifically, we generate Ns samples from a zero-mean Gaussian with
% information matrix A by repeatedly solving the linear system
% A*x_s = b_s, s=1:Ns, 
% with
% b_s = X'*R*mu_rs + B'*P*mu_ps
% and
% mu_rs ~ Gauss(0,inv(R)),   mu_ps ~ Gauss(0,inv(P))
%
% Then
% z \approx 1/Ns*sum_s (B*x_s).^2
%
% ARGUMENTS
% X      [m,n]  matrix or operator, m can be zero
% R      [m,m]  square SPD inverse covariance matrix
%     or [m,1]  vector of diagonal          R = diag(R(:))
%     or [1,1]  scalar multiple of identity R = R*eye(n)
%
% B      [q,n]  matrix or operator
% P      [q,q]  square SPD inverse covariance matrix
%     or [q,1]  vector of diagonal          P = diag(P(:))
%     or [1,1]  scalar multiple of identity P = P*eye(n)
%
% Ns     [1,1]  number of samples           [default 20]
% Ncg    [1,1]  number of LCG iterations    [default 100]
% 
% RESULTS
% z      [q,1]  Unbiased Monte-Carlo estimate of diag(B*inv(A)*B'), 
%               A = X'*R*X + B'*P*B. The estimator has rank-Ns and its relative
%               error drops as 1/sqrt(Ns)
%               It also enforces the constraint that z <= diag(inv(P))
% zu     [n,1]  Unbiased Monte-Carlo estimate of diag(inv(A)), 
%               A = X'*R*X + B'*P*B. The estimator has rank-Ns and its relative
%               error drops as 1/sqrt(Ns)
% ldA    [1,1]  Monte-Carlo estimate of the log determinant of A. It works
%               by comparing the logdet(A) with the logdet(M), where M is
%               the circular preconditioner used in diaginv_sample. It is
%               reliable when the preconditioner is reasonably close to A.
% Q, T          Empty output arguments to obey the interface of diaginv_*.m
%
% The matrices X and B can be given implicitly through their multiplication.
%
% G. Papandreou and A. Yuille, Gaussian sampling by local perturbations, NIPS-10
%
% (c) by George Papandreou, CIVS@UCLA, 2011 October 21

function [z,zu,ldA,Q,T] = diaginv_sample(X,R,B,P,Ns,Ncg)

if nargin<5, Ns = 100; end                                   % number of samples
if nargin<6, Ncg = 100; end

Rsqrt = info_sqrt(R);
Psqrt = info_sqrt(P);

% Find sizes
n = max(size(X,2),size(B,2));
if isempty(X)
  m = 0;
elseif isscalar(X)
  m = max(size(R));
else
  m = size(X,1);
end
if isempty(B)
  error('diaginv_sample: B cannot be empty');
elseif isscalar(B)
  q = max(size(P));
else
  q = size(B,1);
end

z = zeros([q 1]);
Es = zeros(Ns,1);
comp_zu = nargout>1;
comp_ldA = nargout>2;

if comp_zu
  zu = zeros(n,1);
else
  zu = [];
end

% logdet(M) has the same value at each iteration, but parfor complains
ldM = zeros(Ns,1);

% in this loop most of the expensive computation takes place
% if you are using a recent version of Matlab, you can replace "for" by "parfor"
% to take advantage of your multiprocessor infrastructure
for s=1:Ns
  k_s = genks(X,Rsqrt,m) + genks(B,Psqrt,q);
  [x_s,cg_done,cg_res,mvmMfun,prec] = linsolve_lcg(X,R,B,P,k_s,Ncg);
  z = z + (B*x_s).^2;
  if comp_zu
    zu = zu + x_s.^2; 
  end
  
  % Es(s) = 0.5*(x_s'*(k_s-x_s));                                % energy if M=I
  [mx_s d] = mvmMfun(x_s);                          % mx_s=M*x_s, M=F'*diag(d)*F
  ldM(s) = sum(log(d));
  Es(s) = 0.5*(x_s'*(k_s-mx_s));
end

z = 1/Ns*z;
z = cov_constraint(P,z);

if comp_ldA
  ldA = ldM(1) + ...                % logdet(M), with M being the preconditioner
        LSum(Es,1) - log(Ns);           % term comparing the proximity of A to M
  if ~prec, warning(['diaginv_sample: ldA, the estimate of the log ',        ...
                     'determinant of A is probably unreliable since no ',    ...
                     'preconditioning is possible'])
  end
end

if comp_zu
  zu = 1/Ns*zu;
end

if nargout>3, Q=[]; end, if nargout>4, T=[]; end

% J      [m,m]  square SPD inverse covariance matrix
%     or [m,1]  vector of diagonal          J = diag(J(:))
%     or [1,1]  scalar multiple of identity J = J*eye(m)
function Jsqrt = info_sqrt(J)
[m n] = size(J);
if isempty(J)
  Jsqrt = [];
elseif m==1 || n==1, 
  Jsqrt = sqrt(J);
elseif m==n
  Jsqrt = chol(J,'lower');
else
  error('diaginv_sample: wrong J');
end

% generates k_s = F'*C*mu_s = F'*C*C^(-0.5)*n_s = F'*C^(0.5)*n_s,
% with mu_s ~ Gauss(0,inv(C)) and n_s = Gauss(0,I)
function k_s = genks(F,Csqrt,d)
if isempty(F)
  k_s = 0;
  return;
end
n_s = randn([d 1]);
[m n] = size(Csqrt);
if m==1 || n==1, 
  z = Csqrt.*n_s;
elseif m==n
  z = Csqrt*n_s;
else
  error('diaginv_sample: wrong Csqrt');
end
k_s = [F']*z;

% Enforce the constraint that z <= diag(inv(P)), whenever this is cheap to
% do, i.e., J is spherical or diagonal
function z = cov_constraint(J,z)
[m n] = size(J);
if  m==1 && n==1
  z = min(z,1/J);
elseif m==1 || n==1, 
  z = min(z,1./reshape(J,size(z)));
end

% returns log(exp(x1)+exp(x2)+...+exp(xn)) with precautions.
% the equivalent of sum, but replaces addition with LAdd
% useful for properly adding log probs
function ls = LSum(x,dim)

if nargin<2, dim=1; end
LZERO = -1e10;
sz = size(x);
 
sz1 = [sz(1:dim-1), sz(dim+1:end)];
if length(sz1) ~= length(sz)
  sz1 = [sz1(:)', 1];
end

sz2 = [sz(1:dim-1), 1, sz(dim+1:end)];
ls = reshape(LZERO*ones(sz1),sz2);

x_ndims = ndims(x);
switch x_ndims
  case 1
    colons = {':'};
  case 2
    colons = {':' ':'};
  case 3
    colons = {':' ':' ':'};
  case 4
    colons = {':' ':' ':' ':'};
  otherwise
    colons = repmat({':'}, [1, x_ndims]);
end
 
for i = 1:size(x,dim)
  colons{dim} = i;
  ls = LAdd(ls,x(colons{:}));
end

% returns log(exp(x)+exp(y)) with precautions.
% used to properly add log-probabilities
function res = LAdd(x,y)

LZERO     = -1e10;
LSMALL    =  -0.5e10;
minLogExp = -23.0259; % = -log(-LZERO);
 
I = x<y;
temp = x(I);
x(I) = y(I);
y(I) = temp;
diff = y-x;

I1 = diff<minLogExp;
temp =  x(I1);

J = temp<LSMALL;
temp(J) = LZERO;
x(I1) = temp;

I2 = diff>=minLogExp;
z = exp(diff(I2));
z = x(I2)+log(1+z);

res = zeros(size(x));
res(I1) = x(I1);
res(I2) = z;