% Penalised least squares solved by TRUNCATED NEWTON 
% i.e. Newton steps approximated by conjugate gradients
%
%   [u,phi] = plsTN(u,X,y,B,opt,lam,pen,varargin)
%
% Additional options:
%  opt.
%      nit:       maximum number of Newton steps                   [default  15]
%                 we have nMVM/nit MVMs to approximate the Newton dir. by CG
%      nline:     max. bisection steps in Brent's line search      [default 10]
%      exactNewt: flag indicating whether we want to compute the Newton
%                 direction exactly instead of CG                [default false]
%                 this works only if B=1 and X is numeric
%
%   See also PLSSOLVERS.M.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 07

function [u,phi] = plsTN(u,X,y,B,opt,lam,pen,varargin)

% optimisation parameters
if isfield(opt,'nMVM')                           % maximum number of MVMs with A
  nMVM = opt.nMVM;
else
  nMVM = 100;
end
if isfield(opt,'nit')                           % maximum number of Newton steps 
  nit = opt.nit;
else
  nit = 15;
end
ncg = max(fix(nMVM/nit),1);    % maximum number of CG steps for Newton direction
if isfield(opt,'nline')            % max. bisection steps in Brent's line search
  nline = opt.nline;
else
  nline = 10;
end
if isfield(opt,'exactNewt')     % Do we want to compute the Newton dir. exactly?
  exactNewt = opt.exactNewt;
else
  exactNewt = false;
end

% information parameters
if isfield(opt,'output')              % flag saying whether some output is shown
  output = opt.output;
else
  output = false;
end

n = numel(u); 
if output, if exist('fflush','builtin'), fflush(stdout); end, end
for i=1:nit
  s = B*u;  r = X*u-y;                                % sparse coeffs, residuals
  [h,dh,d2h] = feval(pen,s,varargin{:}); g = [X']*r/lam+[B']*dh;
  if exactNewt
    woodbury = isnumeric(B) && numel(B)==1 && isnumeric(X);
    if woodbury, woodbury = woodbury && B==1; end        % Is Woodbury possible?
    if woodbury
      % compute Newton direction by exploiting the Woodbury formula d = -H\g;
      d = linsolve_woodbury( X, 1/lam, d2h, -g );
    else
      % compute Newton direction by full matrices
      d = linsolve_full( X, 1/lam, B, d2h, -g);
    end
  else
    % compute Newton direction by LCG: d = -H\g;
    %  gradient g = X'*(X*u-y)/lam + B'*dh;
    %  Hessian  H = X'*X/lam       + B'*diag(d2h)*B;
    %  A = X'*R*X + B'*diag(pi)*B, R = eye(m)/lam, P = diag(d2h)
    d = linsolve_lcg( X, 1/lam, B, d2h, -g,min(ncg,n+1));
  end

  % line search along Newton direction d, stepsize t
  Phi = @(t,r,X_d,s,B_d,lam,varargin) ...
                  norm(r+t*X_d)^2/lam + 2*sum(feval(pen,s+t*B_d,varargin{:}));
  [t,phi] = brentmin(0,2,nline,1e-8,Phi,r,X*d,s,B*d,lam,varargin{:});  % t=[0,2]

  % update estimate
  u = u+t*d; du = t*norm(d)/max(norm(u),1e-9);
  if output
    fprintf('%5i, phi=%4.4e;  du=%1.3e\r', i, phi, du )
    if exist('fflush','builtin'), fflush(stdout); end
  end
  if du<sqrt(eps), break, end                                        % converged
end
if output, fprintf('\n'); if exist('fflush','builtin'), fflush(stdout); end, end
