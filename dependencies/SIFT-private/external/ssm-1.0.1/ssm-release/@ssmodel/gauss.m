function [model ytilde alpha] = gauss(y, model, alpha, varargin)

%@SSMODEL/GAUSS Calculate Gaussian approximation of model.
%   [model ytilde alpha] = GAUSS(y, model[, alpha0])
%   GAUSS(..., optname1, optvalue1, optname2, optvalue2, ...)
%       y is the observation sequence.
%       model is the state space model to use.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3 || isempty(alpha), alpha = zeros(size(model.T, 1), size(y, 2));
else checkstate(model, alpha, 'alpha0'); end

opt     = setopt([], varargin{:});

checkdata(model, y, 'y');

if opt.usec
    %% Gaussian approximation %%
    [model.H model.Q ytilde alpha converged c] = gauss_int(size(y, 2), y, [], [], [], isa(model.H, 'ssdist'), isa(model.Q, 'ssdist'), model.H, model.Q, ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), getmat_c(model.Z), getmat_c(model.T), getmat_c(model.R), getmat_c(model.c), getmat_c(model.a1), getmat_c(model.P1), alpha, opt.tol, opt.maxiter, true, opt.inv);
else
    %% Data preprocessing %%
    mis = isnan(y);

    %% Gaussian approximation %%
    [model.H model.Q ytilde alpha converged c] = gauss_int(size(y, 2), y, mis, any(mis, 1), all(mis, 1), isa(model.H, 'ssdist'), isa(model.Q, 'ssdist'), model.H, model.Q, ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.c), model.a1.mat, model.P1.mat, alpha, opt.tol, opt.maxiter, false, opt.inv);
end

%% Output display %%
if opt.disp >= 1 && ~converged, fprintf(1, 'Warning: Gaussian approximation didn''t converge.\n'); end
if opt.disp >= 2
    if converged, fprintf(1, 'Gaussian Approximation Result: model converged in %d iterations.\n', c);
    else fprintf(1, 'Gaussian Approximation Result: model didn''t converge.\n'); end
end

