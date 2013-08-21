% lbfgs - L-BFGS algorithm
%
% Syntax:
%  [xx, status] = lbfgs(fun, xx, ll, uu, <opt>)
%
% Input:
%  fun     - objective function
%  xx      - Initial point for optimization
%  ll      - lower bound on xx
%  uu      - upper bound on xx
%  Ac      - inequality constraint:
%  bc      -     Ac*xx<=bc
%  opt     - Struct or property/value list of optional properties:
%   .m          - size of limited memory
%   .epsg       - gradient tolerance
%   .maxiter    - maximum number of iterations
%   .display    - display level
% 
% Output:
%  xx      - Final point of optimization
%  status  - Various numbers
%
% Copyright(c) 2009 Ryota Tomioka
% This software is distributed under the MIT license. See license.txt

function [xx, status] = lbfgs(fun, xx, ll, uu, Ac, bc, info, opt, varargin)

opt=set_defaults(opt, 'm', 6,...
                      'ftol', 1e-5, ...
                      'maxiter', 0,...
                      'max_linesearch', 50,...
                      'display', 0,...
                      'epsg', 1e-5,...
                      'epsginfo', 1);

if ischar(fun)
  fun = {fun};
end

nn = size(xx,1);

t0 = cputime;

% Limited memory
lm = repmat(struct('s',zeros(nn,1),'y',zeros(nn,1),'ys',0,'alpha',0),[1, opt.m]);

[fval,gg,info]=fun(xx, info);


% The initial step is gradient
dd = -gg;

kk = 1;
stp = 1/norm(dd);

bResetLBFGS = 0;
ixend = 1;
bound = 0;
while 1
  fp = fval;
  xxp = xx;
  ggp = gg;

  % Perform line search
  [ret, xx,fval,gg,info,stp]=...
      linesearch_backtracking(fun, xx, ll, uu, Ac, bc, fval, gg, dd, stp, info, opt, varargin{:});
  
  if ret<0
    fprintf('ginfo=%g\n',info.ginfo);
    break;
  end

  % Progress report
  gnorm = norm(gg);
  if opt.display>1
    fprintf('[%d] xx=[%g %g...] fval=%g gnorm=%g step=%g\n',kk,xx(1),xx(2),fval,gnorm,stp);
  end

  if info.ginfo<opt.epsginfo % || gnorm<opt.epsg
    if opt.display>1
      fprintf('Optimization success! ginfo=%g\n',info.ginfo);
    end
    ret=0;
    break;
  end
  
  if kk==opt.maxiter
    if opt.display>0
      fprintf('Maximum #iterations=%d reached.\n', kk);
    end
    ret = -3;
    break;
  end

  % L-BFGS update
  if opt.m>0
    lm(ixend).s = xx-xxp;
    lm(ixend).y = gg-ggp;
    ys = lm(ixend).y'*lm(ixend).s; yy = sum(lm(ixend).y.^2);
    lm(ixend).ys  = ys;
  else
    ys = 1; yy = 1;
  end
  
  bound = min(bound+1, opt.m);
  ixend = (opt.m>0)*(mod(ixend, opt.m)+1);

  % Initially set the negative gradient as descent direction
  dd = -gg;
  
  jj = ixend;
  for ii=1:bound
    jj = mod(jj + opt.m -2, opt.m)+1;
    lm(jj).alpha = lm(jj).s'*dd/lm(jj).ys;
    dd = dd -lm(jj).alpha*lm(jj).y;
  end

  dd = dd *(ys/yy);
  
  for ii=1:bound
    beta = lm(jj).y'*dd/lm(jj).ys;
    dd = dd + (lm(jj).alpha-beta)*lm(jj).s;
    jj = mod(jj,opt.m)+1;
  end

  stp = 1.0;
  
  kk = kk + 1;
end

status=struct('ret', ret,...
              'kk', kk,...
              'fval', fval,...
              'gg', gg,...
              'time', cputime-t0,...
              'info', info,...
              'opt', opt);


function [ret, xx, fval, gg, info, step]...
    =linesearch_backtracking(fun, xx, ll, uu, Ac, bc, fval, gg, dd, step, info, opt, varargin)

floss=0;
gloss=zeros(size(gg));

dginit=gg'*dd;

if dginit>=0
  if opt.display>0
    fprintf('dg=%g is not a descending direction!\n', dginit);
  end
  step = 0;
  ret = -1;
  return;
end

Ip=find(dd>0);
In=find(dd<0);
step=min([step, 0.999*min((xx(In)-ll(In))./(-dd(In))), 0.999*min((uu(Ip)-xx(Ip))./dd(Ip))]);


xx0 = xx;
f0  = fval;
gg0 = gg;
cc = 0;

if opt.display>2
  fprintf('finit=%.20f\n',f0);
end

while cc<opt.max_linesearch
  ftest = f0  + opt.ftol*step*dginit;
  xx    = xx0 + step*dd;

  if ~isempty(Ac)
    bineq = all(Ac*xx<=bc);
  else
    bineq = true;
  end

  if bineq && all(xx>=ll) && all(xx<=uu)
    [fval, gg, info]=fun(xx, info);
    
    if fval<=ftest
      break;
    end
  else
    fval = inf;
  end
  if opt.display>2
    fprintf('[%d] step=%g fval=%.20f > ftest=%.20f\n', cc, step, fval, ftest);
  end
  
  step = step/2;
  cc = cc+1;
end

if cc==opt.max_linesearch
  if opt.display>0
    fprintf('Maximum linesearch=%d reached\n', cc);
  end
  xx   = xx0;
  [fval, gg, info]=fun(xx, info);
  step = 0;
  ret = -2;
  return;
end


ret = 0;
