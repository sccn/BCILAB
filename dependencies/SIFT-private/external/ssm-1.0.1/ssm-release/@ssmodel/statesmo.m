function [alphahat V r N] = statesmo(y, model, varargin)

%@SSMODEL/STATESMO State smoothing with the linear Gaussian model.
%   [alphahat V r N] = STATESMO(y, model)
%   STATESMO(..., optname1, optvalue1, optname2, optvalue2, ...)
%       y is the observation sequence.
%       model is the state space model to use.
%       alphahat is the smoothed state (expected state given all observation).
%       V is the smoothed state variance (state variance given all observation).

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

opt     = setopt([], varargin{:});

checkdata(model, y, 'y');

if opt.usec
    if nargout > 1
        [alphahat V r N] = statesmo_int_c(y, getmat_c(model.H), ~issta(model.H), getmat_c(model.Z), ~issta(model.Z), ...
            getmat_c(model.T), ~issta(model.T), getmat_c(model.R), ~issta(model.R), getmat_c(model.Q), ~issta(model.Q), ...
            getmat_c(model.c), ~issta(model.c), getmat_c(model.a1), getmat_c(model.P1), opt.tol, opt.tol, opt.inv, true);

        %% Result postprocessing %%
        if ndims(y) > 2, r(:,:, 1) = []; else r(:, 1) = []; end
        if ndims(N) > 2, N(:,:, 1) = []; else N(:, 1) = []; end
    else
        alphahat = fastsmo_int_c(y, getmat_c(model.H), ~issta(model.H), getmat_c(model.Z), ~issta(model.Z), ...
            getmat_c(model.T), ~issta(model.T), getmat_c(model.R), ~issta(model.R), getmat_c(model.Q), ~issta(model.Q), ...
            getmat_c(model.c), ~issta(model.c), getmat_c(model.a1), getmat_c(model.P1), opt.tol, opt.tol, opt.inv, true);
    end
else
    %% Data preprocessing %%
    n   = size(y, 2);
    mis = isnan(y);

    if ndims(y) > 2
        if any(mis(:)), warning('ssm:ssmodel:statesmo:BatchMissing', 'Batch operations with missing data not supported for MATLAB code, set ''usec'' to true.'); end

        %% Data preprocessing %%
        N   = size(y, 3);
        y   = permute(y, [1 3 2]);

        if nargout > 1
            %% Batch smoother %%
            [alphahat temp temp r] = batchsmo_int(n, N, y, ...
                ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), ...
                getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), ...
                model.a1.mat, model.P1.mat, opt.tol);

            %% Get V and N %%
            [temp V temp N] = statesmo_int(1, n, y, mis, any(mis, 1), all(mis, 1), ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), model.a1.mat, model.P1.mat, opt.tol);

            %% Result postprocessing %%
            alphahat    = permute(alphahat, [1 3 2]);
            r           = permute(cat(3, r{:}), [1 3 2]);
            N           = cat(3, N{:});
            if size(V, 1) == 1
                V   = permute(V, [2 3 1]);
                N   = permute(N, [2 3 1]);
            end
        else
            %% Batch smoother %%
            alphahat = batchsmo_int(n, N, y, ...
                ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), ...
                getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), ...
                model.a1.mat, model.P1.mat, opt.tol);

            %% Result postprocessing %%
            alphahat = permute(alphahat, [1 3 2]);
        end
    else
        %% State smoothing %%
        if nargout > 1
            [alphahat V r N] = statesmo_int(1, n, y, mis, any(mis, 1), all(mis, 1), ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), model.a1.mat, model.P1.mat, opt.tol);

            %% Result postprocessing %%
            r   = [r{:}];
            N   = cat(3, N{:});
            if size(V, 1) == 1
                V   = permute(V, [2 3 1]);
                N   = permute(N, [2 3 1]);
            end
        else alphahat = faststatesmo_int(n, y, mis, any(mis, 1), all(mis, 1), ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), model.a1.mat, model.P1.mat, opt.tol);
        end
    end
end

