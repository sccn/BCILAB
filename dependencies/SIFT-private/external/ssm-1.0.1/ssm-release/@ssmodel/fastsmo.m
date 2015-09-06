function [alphahat epshat etahat] = fastsmo(y, model, varargin)

%@SSMODEL/FASTSMO Fast smoothing with the linear Gaussian model.
%   [alphahat epshat etahat] = FASTSMO(y, model)
%   FASTSMO(..., optname1, optvalue1, optname2, optvalue2, ...)
%       y is the observation sequence.
%       model is the state space model to use.
%       alphahat is the smoothed state (expected state given all observation).
%       epshat is the smoothed observation disturbance.
%       etahat is the smoothed state disturbance.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

opt     = setopt([], varargin{:});

checkdata(model, y, 'y');

if opt.usec
    [alphahat epshat etahat] = fastsmo_int_c(y, getmat_c(model.H), ~issta(model.H), getmat_c(model.Z), ~issta(model.Z), ...
        getmat_c(model.T), ~issta(model.T), getmat_c(model.R), ~issta(model.R), getmat_c(model.Q), ~issta(model.Q), ...
        getmat_c(model.c), ~issta(model.c), getmat_c(model.a1), getmat_c(model.P1), opt.tol, opt.tol, opt.inv, true);
else
    %% Data preprocessing %%
    mis = isnan(y);

    if ndims(y) > 2
        if any(mis(:)), warning('ssm:ssmodel:fastsmo:BatchMissing', 'Batch operations with missing data not supported for MATLAB code, set ''usec'' to true.'); end

        %% Data preprocessing %%
        n   = size(y, 2);
        N   = size(y, 3);
        y   = permute(y, [1 3 2]);

        %% Batch smoother %%
        [alphahat epshat etahat] = batchsmo_int(n, N, y, ...
            ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), ...
            getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), ...
            model.a1.mat, model.P1.mat, opt.tol);

        %% Result postprocessing %%
        alphahat = permute(alphahat, [1 3 2]);
        epshat   = permute(epshat, [1 3 2]);
        etahat   = permute(etahat, [1 3 2]);
    else
        %% Fast smoother %%
        [alphahat epshat etahat] = fastsmo_int(size(y, 2), y, mis, any(mis, 1), all(mis, 1), ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), model.a1.mat, model.P1.mat, opt.tol);
    end
end

