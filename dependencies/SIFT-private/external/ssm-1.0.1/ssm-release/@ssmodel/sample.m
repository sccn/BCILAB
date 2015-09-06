function [y alpha eps eta] = sample(model, n, N, varargin)

%@SSMODEL/SAMPLE Generate unconditional samples from linear Gaussian model.
%   [y alpha eps eta] = SAMPLE(model, n[, N])
%   SAMPLE(..., optname1, optvalue1, optname2, optvalue2, ...)
%       model is the linear Gaussian model to use.
%       n is the time duration to generate.
%       N is the number of samples to generate.
%       y is the p*n*N observation samples generated.
%       alpha, eps, eta the corresponding .*n*N state and disturbances.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, error('ssm:ssmodel:sample:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(n) || ~isnumeric(n), error('ssm:ssmodel:sample:InputError', 'n must be a scalar.'); end
if nargin < 3, N = 1; elseif ~isscalar(N) || ~isnumeric(N), error('ssm:ssmodel:sample:InputError', 'N must be a scalar.'); end

opt     = setopt([], varargin{:});

Znl     = isa(model.Z, 'ssfunc');
Tnl     = isa(model.T, 'ssfunc');

if opt.usec && ~Znl && ~Tnl
    [y alpha eps eta] = sample_int_c(n, N, getmat_c(model.H), ~issta(model.H), getmat_c(model.Z), ~issta(model.Z), ...
        getmat_c(model.T), ~issta(model.T), getmat_c(model.R), ~issta(model.R), getmat_c(model.Q), ~issta(model.Q), ...
        getmat_c(model.c), ~issta(model.c), getmat_c(model.a1), getmat_c(model.P1), true);
else
    %% Unconditional sampling %%
    [y alpha eps eta] = sample_int(n, N, model.p, size(model.T, 1), size(model.R, 2), ...
        isa(model.Z, 'ssfunc'), isa(model.T, 'ssfunc'), model.Z, model.T, ...
        ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), ...
        getmat(model.H), getmat(model.R), getmat(model.Q), getmat(model.c), model.a1.mat, model.P1.mat);

    %% Result postprocessing %%
    y       = permute(y, [1 3 2]);
    alpha   = permute(alpha, [1 3 2]);
    eps     = permute(eps, [1 3 2]);
    eta     = permute(eta, [1 3 2]);
end

