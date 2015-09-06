function [a P v F] = kalman(y, model, varargin)

%@SSMODEL/KALMAN Kalman filtering with the linear Gaussian model.
%   [a P v F] = KALMAN(y, model) kalman filters y using model.
%   KALMAN(..., optname1, optvalue1, optname2, optvalue2, ...)
%       y is the observation sequence.
%       model is the state space model to use.
%       a is the expected state given previous observation.
%       P is the state variance given previous observation.
%       v is the one-step-ahead forecast error.
%       F is the one-step-ahead forecast error variance.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

opt     = setopt([], varargin{:});

checkdata(model, y, 'y');

if opt.usec
    [a P v F] = kalman_int_c(y, getmat_c(model.H), ~issta(model.H), getmat_c(model.Z), ~issta(model.Z), ...
        getmat_c(model.T), ~issta(model.T), getmat_c(model.R), ~issta(model.R), getmat_c(model.Q), ~issta(model.Q), ...
        getmat_c(model.c), ~issta(model.c), getmat_c(model.a1), getmat_c(model.P1), opt.hideinit, opt.tol, opt.tol, opt.inv, true);
else
    %% Data preprocessing %%
    p       = size(y, 1);
    n       = size(y, 2);
    mis     = isnan(y);
    anymis  = any(mis, 1);
    allmis  = all(mis, 1);

    if ndims(y) > 2
        if any(mis(:)), warning('ssm:ssmodel:kalman:BatchMissing', 'Batch operations with missing data not supported for MATLAB code, set ''usec'' to true.'); end

        N   = size(y, 3);
        y   = permute(y, [1 3 2]);

        %% Kalman filter %%
        [a P d vc invFc]    = batchkalman_int(1, n, N, y, ...
            ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), ...
            getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), ...
            getmat(model.a1), getmat(model.P1), opt.tol);

        %% Result postprocessing %%
        if opt.hideinit
            a(:,:, 1:d) = NaN;
            P(:,:, 1:d) = NaN;
        end
        a   = permute(a, [1 3 2]);
        if size(P, 1) == 1, P = permute(P, [2 3 1]); end % simplify P's representation if m = 1
        if nargout > 2
            v       = zeros(p, N, n);
            F       = zeros(p, p, n);
            if opt.hideinit
                v(:,:, 1:d) = NaN;
                F(:,:, 1:d) = NaN;
            end
            for t = d+1 : n
                if anymis(t) && ~allmis(t)
                    v(~mis(:, t), :, t)             = vc{t};
                    F(~mis(:, t), ~mis(:, t), t)    = inv(invFc{t});
                else
                    v(:, :, t)  = vc{t};
                    F(:, :, t)  = inv(invFc{t});
                end
            end
            if p == 1, F = permute(F, [2 3 1]); end
            v   = permute(v, [1 3 2]);
        end
    else
        %% Kalman filter %%
        [a P d vc invFc]    = kalman_int(1, n, y, mis, anymis, allmis, ...
            ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), ...
            getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), ...
            getmat(model.a1), getmat(model.P1), opt.tol);

        %% Result postprocessing %%
        if opt.hideinit
            a(:, 1:d)   = NaN;
            P(:,:, 1:d) = NaN;
        end
        if size(P, 1) == 1, P = permute(P, [2 3 1]); end % simplify P's representation if m = 1
        if nargout > 2
            v       = zeros(p, n);
            F       = zeros(p, p, n);
            if opt.hideinit
                v(:, 1:d)   = NaN;
                F(:,:, 1:d) = NaN;
            end
            for t = d+1 : n
                if anymis(t) && ~allmis(t)
                    v(~mis(:, t), t)                = vc{t};
                    F(~mis(:, t), ~mis(:, t), t)    = inv(invFc{t});
                else
                    v(:, t)     = vc{t};
                    F(:, :, t)  = inv(invFc{t});
                end
            end
            if p == 1, F = permute(F, [2 3 1]); end
        end
    end
end

