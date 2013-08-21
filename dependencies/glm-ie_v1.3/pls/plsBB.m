% Penalised least squares by Barzilai/Borwein stepsize gradient descent [1]
%
%   [u,phi] = plsBB(u,X,y,B,opt,lam,pen,varargin)
%  
%   See also PLSSOLVERS.M.
%
%   The optimiser is purely gradient based; it does not use the objective 
%   function values for line searches etc. 
%   Further, the solver is not guaranteed to decrease the objective in every
%   step. The solver has very good empirical behavior for high-dimensional
%   optimisation problems such as l1-regularised least squares.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 January 19
%
% [1] Jonathan Barzilai and Jonathan M. Borwein, Two-Point Step Size Gradient
% Methods, IMA Journal of Numerical Analysis (1988) 8, 141-148

function [u,phi] = plsBB(u,X,y,B,opt,lam,pen,varargin)

% optimisation parameters
if isfield(opt,'eps_grad')                 % gradient norm convergence threshold
  eps_grad = opt.eps_grad;
else
  eps_grad = 1e-10;
end
if isfield(opt,'nMVM')                                 % maximal number of steps
  nMVM = opt.nMVM;
else
  nMVM = 100;
end
if isfield(opt,'plsBBstepsizetype')                          % type of step size
  plsBBstepsizetype = opt.plsBBstepsizetype;
else
  plsBBstepsizetype = 1;
end

% information parameters
if isfield(opt,'output')              % flag saying whether some output is shown
  output = opt.output;
else
  output = false;
end

if output, if exist('fflush','builtin'), fflush(stdout); end, end
for n=1:nMVM
  [f,g] = feval('phi',u,X,y,B,lam,pen,varargin{:});
  if n==1
    uold = randn(size(u))/1e2;                 % draw a virtual initial location
    [fold,gold] = feval('phi',uold,X,y,B,lam,pen,varargin{:});
    if fold<f % swap
      tmp = gold; gold = g; g = tmp;
      tmp = uold; uold = u; u = tmp;
                  fold = f;
    end
  end
  dg = g-gold; du = u-uold;
  if plsBBstepsizetype==1
    a = (du'*dg)/(dg'*dg);
  else
    a = (du'*du)/(du'*dg);
  end
  uold = u; gold = g; u = u - a*g;       % update cache and variable of interest
  if output
    fprintf('%5i, phi=%4.4e;  du=%1.3e\r', n, f, norm(du)/max(norm(u),1e-9) )
    if exist('fflush','builtin'), fflush(stdout); end
  end
  if norm(du(:)) < eps_grad, break, end                     % stopping criterion
end
if output, fprintf('\n'); if exist('fflush','builtin'), fflush(stdout); end, end

phi = feval('phi',u,X,y,B,lam,pen,varargin{:});