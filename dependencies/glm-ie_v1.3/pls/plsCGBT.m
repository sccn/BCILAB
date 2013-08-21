% Penalised least squares solved by CONJUGATE GRADIENTS with BACKTRACKING line 
% search using the Armijo rule. The algorithm is inspired by the code fnlCg.m 
% from Michael Lustig http://www.stanford.edu/~mlustig/sparseMRI_v0.2.tar.gz.
%
%   [u,phi] = plsTN(u,X,y,B,opt,lam,pen,varargin)
%
% Additional options:
%  opt.
%      nline      : maximum number of line searche steps        [default  20]
%      alpha      : constant for Armijo condition               [default 0.01]
%      beta       : stepsize update factor                      [default 0.6]
%      cmax       : max. stepsize (gradient scaling factor)     [default 1]
%      eps_grad   : gradient norm convergence threshold         [default 1e-10]
%      nrestart   : number of times CG is restarted             [default 0]
%
%   See also PLSSOLVERS.M.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 07

function [u,phi] = plsCGBT(u,X,y,B,opt,lam,pen,varargin)

% line search parameters
if isfield(opt,'nline')                   % maximum number of line searche steps
  nline = opt.nline;
else
  nline = 20;
end
if isfield(opt,'alpha')                          % constant for Armijo condition
  al = opt.alpha;
else
  al = 0.01;
end
if isfield(opt,'beta')                                  % stepsize update factor
  be = opt.beta;
else
  be = 0.6;
end
if isfield(opt,'cmax')                 % max. stepsize (gradient scaling factor)
  cmax = opt.cmax;
else
  cmax = 1.0; opt.cmax = 1.0;
end

% optimisation parameters
if isfield(opt,'eps_grad')                 % gradient norm convergence threshold
  eps_grad = opt.eps_grad;
else
  eps_grad = 1e-10;
end
if isfield(opt,'nMVM')              % maximal number of conjugate gradient steps
  nMVM = opt.nMVM;
else
  nMVM = 100;
end
if isfield(opt,'nrestart')                               % number of CG restarts
  nrestart = opt.nrestart;
  idrestart = fix(linspace(0,nMVM,nrestart+2)); 
  idrestart = idrestart(1:end-1);                              % when to restart
else
  idrestart = 0;
end

% information parameters
if isfield(opt,'output')              % flag saying whether some output is shown
  output = opt.output;
else
  output = false;
end

varargin = {lam,pen,varargin{:}};
if output, if exist('fflush','builtin'), fflush(stdout); end, end
for n=0:nMVM
  if any(n==idrestart)                                                 % restart
    [g0,Xu,Bu] = gradPhi(u,X,y,B,varargin{:}); du = -g0;      % initial gradient
    cmax = opt.cmax;
  end
  Xdu = X*du; Bdu = B*du;                         % pre-calculate MVMs for speed
	c=0;    f0 = Phi(Xu,Xdu,y,Bu,Bdu,c,varargin{:});           % objective for c=0
  c=cmax; f1 = Phi(Xu,Xdu,y,Bu,Bdu,c,varargin{:});              % and for c=cmax 
  for nl = 0:nline-1                                  % backtracking line-search
    if f1 <= f0-c*al*abs(g0(:)'*du(:)), break, end      % check Armijo condition
    c = c*be;  f1 = Phi(Xu,Xdu,y,Bu,Bdu,c,varargin{:});
  end
  % control the number of line searches by adapting the initial step search
	if nl>2, cmax = cmax*be; end, if nl<1, cmax = cmax/be; end
  if nl==nline, warning(1,'Max line search; operator bug?'); return, end
	u  = u + c*du;                                                        % update
	[g1,Xu,Bu] = gradPhi(u,X,y,B,varargin{:});         % Fletcher-Reeves CG update
	bk = g1(:)'*g1(:)/( g0(:)'*g0(:)+eps );  g0 = g1;  du = bk*du-g1;
  if output
    fprintf('%5i, phi=%4.4e;  du=%1.3e\r', n, f1, c*norm(du)/max(norm(u),1e-9) )
    if exist('fflush','builtin'), fflush(stdout); end
  end
  if norm(du(:)) < eps_grad, break, end                     % stopping criterion
end
phi = f1;
if output, fprintf('\n'); if exist('fflush','builtin'), fflush(stdout); end, end


function phi = Phi(Xu,Xdu,y,Bu,Bdu,c,lam,pen,varargin)   % objective Phi(u+c*du)
  s = Bu+c*Bdu;  r = Xu+c*Xdu - y;
  h = feval(pen,s,varargin{:});
  phi = (r'*r)/lam + 2*sum(h);


function [g,Xu,Bu] = gradPhi(u,X,y,B,lam,pen,varargin)   % gradient of objective
  Bu = B*u; Xu = X*u; r = Xu-y;
  [h,dh] = feval(pen,Bu,varargin{:});
  g = 2*([X']*r + lam*([B']*dh));