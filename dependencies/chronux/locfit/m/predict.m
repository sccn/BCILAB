function [y, se] = predict(varargin)

% Interpolate a fit produced by locfit().
%
% predict(fit)    produces the fitted values at locfit's selected points.
% predict(fit,x)  interpolates the fits to points specified by x.
%
% Input arguments:
%   fit   The locfit() fit.
%   x     Points to interpolate at. May be a matrix with d columns,
%         or cell with d components (each a vector). In the former
%         case, a fitted value is computed for each row of x.
%         In the latter, the components of x are interpreted as
%         grid margins.
%         Can also specify 'data' (evaluate at data points);
%         or 'fitp' (extract the fitted points).
%  'band',value
%         Type of standard errors to compute. Default is 'band','n', for none.
%         Other choices are 'band','g' (use a global s to estimate the resiudal
%         standard deviation, so standard errors are s*||l(x)||);
%         'band','l' (use a local s(x), so std. errors are s(x)*||l(x)||);
%         'band','p' (prediction errors, so s*sqrt(1+||l(x)||^2).
%  'direct'
%         Compute the local fit directly (rather than using local
%         regression, at each point specified by the x argument.
%  'kappa',vector
%         Vector of constants for simultaneous confidence bands,
%         computed by the kappa0() function.
%  'level',value
%         Coverage probability for confidence intervals and bands.
%         Default is 0.95.
%
%  Output is a vector of fitted values (if 'band','n'), or a cell
%  with fitted value, standard error vectors, and matrix of lower
%  and upper confidence limits.
%
%  Note that for local likelihood fits, back-transformation is
%  not performed, so that (e.g.) for Poisson regression with the
%  log-link, the output estimates the log-mean, and its standard errors.
%  Likewise, for density estimation, the output is log(density).
%
%  Author: Catherine Loader.

if (nargin<1)
    error('predict requires fit argument');
end;

fit = varargin{1};

if (nargin==1) x = 'fitp'; else x = varargin{2}; end;

band = 'n';
what = 'coef';
rest = 'none';
dir  = 0;
level = 0.95;
d = size(fit.data.x,2);
kap = [zeros(1,d) 1];

na = 3;
while na <= nargin
  inc = 0;
  if strcmp(varargin{na},'band')
    band = varargin{na+1};
    inc = 2;
  end;
  if strcmp(varargin{na},'what')
    what = varargin{na+1};
    inc = 2;
  end;
  if strcmp(varargin{na},'restyp')
    rest = varargin{na+1};
    inc = 2;
  end;
  if strcmp(varargin{na},'direct')
    dir = 1;
    inc = 1;
  end;
  if strcmp(varargin{na},'kappa')
    kap = varargin{na+1};
    inc = 2;
  end;
  if strcmp(varargin{na},'level')
    level = varargin{na+1};
    inc = 2;
  end;
  if (inc == 0)
    disp(varargin{na});
    error('Unknown argument');
  end;
  na = na+inc;
end;

[y se cb] = mexpp(x,fit,band,what,rest,dir,kap,level);
if (band=='n')
    y = y;
else
    y = {y se cb};
end;

return;
