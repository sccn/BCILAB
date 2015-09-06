function [model ytilde alpha] = linear(y, model, alpha, varargin)

%@SSMODEL/LINEAR Calculate linear approximation of model.
%   [model ytilde alpha] = LINEAR(y, model[, alpha0])
%   LINEAR(..., optname1, optvalue1, optname2, optvalue2, ...)
%       y is the observation sequence.
%       model is the state space model to use.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3 || isempty(alpha), alpha = zeros(size(model.T, 1), size(y, 2));
else checkstate(model, alpha, 'alpha0'); end

opt     = setopt([], varargin{:});

checkdata(model, y, 'y');

if opt.usec
    %% Linear approximation %%
    [model.Z model.T model.c ytilde alpha converged c] = linear_int(size(y, 2), y, [], [], [], isa(model.Z, 'ssfunc'), isa(model.T, 'ssfunc'), model.Z, model.T, model.c, ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), getmat_c(model.H), getmat_c(model.R), getmat_c(model.Q), getmat_c(model.a1), getmat_c(model.P1), alpha, opt.tol, opt.maxiter, true, opt.inv);
else
    %% Data preprocessing %%
    mis = isnan(y);

    %% Linear approximation %%
    [model.Z model.T model.c ytilde alpha converged c] = linear_int(size(y, 2), y, mis, any(mis, 1), all(mis, 1), isa(model.Z, 'ssfunc'), isa(model.T, 'ssfunc'), model.Z, model.T, model.c, ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), getmat(model.H), getmat(model.R), getmat(model.Q), model.a1.mat, model.P1.mat, alpha, opt.tol, opt.maxiter, false, opt.inv);
end

%% Output display %%
if opt.disp >= 1 && ~converged, fprintf(1, 'Warning: Linear approximation didn''t converge.\n'); end
if opt.disp >= 2
    if converged, fprintf(1, 'Linear Approximation Result: model converged in %d iterations.\n', c);
    else fprintf(1, 'Linear Approximation Result: model didn''t converge.\n'); end
end

