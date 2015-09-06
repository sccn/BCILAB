function [epshat etahat epsvarhat etavarhat Result_r Result_N] = disturbsmo(y, model, varargin)

%@SSMODEL/DISTURBSMO Disturbance smoothing with the linear Gaussian model.
%   [epshat etahat epsvarhat etavarhat r N] = DISTURBSMO(y, model)
%   DISTURBSMO(..., optname1, optvalue1, optname2, optvalue2, ...)
%       y is the observation sequence.
%       model is the state space model to use.
%       epshat is the smoothed observation disturbance.
%       etahat is the smoothed state disturbance.
%       epsvarhat is the smoothed observation disturbance variance.
%       etavarhat is the smoothed state disturbance variance.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

opt     = setopt([], varargin{:});

checkdata(model, y, 'y');

%% Get the current matrix values %%
Hdyn    = ~issta(model.H);
Zdyn    = ~issta(model.Z);
Tdyn    = ~issta(model.T);
Qdyn    = ~issta(model.Q);

if opt.usec
    if nargout > 2
        [epshat etahat epsvarhat etavarhat Result_r Result_N] = disturbsmo_int_c(y, getmat_c(model.H), Hdyn, ...
            getmat_c(model.Z), Zdyn, getmat_c(model.T), Tdyn, getmat_c(model.R), ~issta(model.R), getmat_c(model.Q), Qdyn, ...
            getmat_c(model.c), ~issta(model.c), getmat_c(model.a1), getmat_c(model.P1), opt.tol, opt.tol, opt.inv, true);
        if ndims(y) > 2, Result_r(:,:, 1) = []; else Result_r(:, 1) = []; end
        if ndims(Result_N) > 2, Result_N(:,:, 1) = []; else Result_N(:, 1) = []; end
    else
        [epshat etahat] = fastsmo_int_c(y, getmat_c(model.H), Hdyn, getmat_c(model.Z), Zdyn, getmat_c(model.T), Tdyn, ...
            getmat_c(model.R), ~issta(model.R), getmat_c(model.Q), Qdyn, getmat_c(model.c), ~issta(model.c), ...
            getmat_c(model.a1), getmat_c(model.P1), opt.tol, opt.tol, opt.inv, true);
    end
else
    Hmat    = getmat(model.H);
    Zmat    = getmat(model.Z);
    Tmat    = getmat(model.T);
    Qmat    = getmat(model.Q);
    a1      = model.a1.mat;
    if ~Hdyn, H = Hmat; end
    if ~Zdyn, Z = Zmat; end
    if ~Tdyn, T = Tmat; end
    if ~Qdyn, Q = Qmat; end

    %% Data preprocessing %%
    n   = size(y, 2);
    mis = isnan(y);
    anymis = any(mis, 1);
    allmis = all(mis, 1);

    if ndims(y) > 2
        if any(mis(:)), warning('ssm:ssmodel:disturbsmo:BatchMissing', 'Batch operations with missing data not supported for MATLAB code, set ''usec'' to true.'); end

        %% Data preprocessing %%
        N   = size(y, 3);
        y   = permute(y, [1 3 2]);

        %% Batch smoother %%
        [temp epshat etahat Result_r] = batchsmo_int(n, N, y, ...
            ~issta(model.H), ~issta(model.Z), ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), ...
            getmat(model.H), getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), ...
            model.a1.mat, model.P1.mat, opt.tol);

        %% Result postprocessing %%
        epshat   = permute(epshat, [1 3 2]);
        etahat   = permute(etahat, [1 3 2]);

        if nargout > 2
            %% Kalman filter %%
            [d Fns v invF K L RQ QRt] = kalman_int(3, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, ~issta(model.R), Qdyn, ~issta(model.c), Hmat, Zmat, Tmat, getmat(model.R), Qmat, getmat(model.c), a1, model.P1.mat, opt.tol);

            %% Disturbance smoothing backwards recursion %%
            p   = size(y, 1);
            m   = size(a1, 1);
            rr  = size(model.Q, 1);
            N   = zeros(m, m);
            epsvarhat   = zeros(p, p, n);
            etavarhat   = zeros(rr, rr, n);
            Result_N    = cell(1, n);
            for t = n : -1 : 1
                Result_N{t}     = N;
                if Hdyn, H = Hmat{t}; end
                if Zdyn, Z = Zmat{t}; end
                if Tdyn, T = Tmat{t}; end
                if Qdyn, Q = Qmat{t}; end
                etavarhat(:,:, t) = Q - QRt{t}*N*RQ{t};
                if allmis(t)
                    %% Disturbance smoothing when all observations are missing %%
                    epsvarhat(:,:, t) = H; %%%% TODO: What is epsvarhat when all missing?
                    N = T'*N*T;
                else
                    if anymis(t), Z(mis(:, t),:)=[]; H(mis(:, t),:)=[]; H(:,mis(:, t))=[]; end
                    if t > d || ~Fns(t)
                        %% Normal disturbance smoothing or when F_inf is zero %%
                        epsvarhat(~mis(:, t), ~mis(:, t), t) = H - H*(invF{t} + K{t}'*N*K{t})*H;
                        M = Z'*invF{t};
                        N = M*Z + L{t}'*N*L{t};
                    else
                        %% Exact initial disturbance smoothing when F_inf is nonsingular %%
                        epsvarhat(~mis(:, t), ~mis(:, t), t) = H - H*K{t}'*N*K{t}*H;
                        N = L{t}'*N*L{t};
                    end
                    if anymis(t), if ~Zdyn, Z = Zmat; end, if ~Hdyn, H = Hmat; end, end
                end
            end
        end
    else
        if nargout > 2
            %% Kalman filter %%
            [d Fns v invF K L RQ QRt] = kalman_int(3, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, ~issta(model.R), Qdyn, ~issta(model.c), Hmat, Zmat, Tmat, getmat(model.R), Qmat, getmat(model.c), a1, model.P1.mat, opt.tol);

            %% Disturbance smoothing backwards recursion %%
            p   = size(y, 1);
            m   = size(a1, 1);
            rr  = size(model.Q, 1);
            r   = zeros(m, 1);
            N   = zeros(m, m);
            epshat      = zeros(p, n);
            epsvarhat   = zeros(p, p, n);
            etahat      = zeros(rr, n);
            etavarhat   = zeros(rr, rr, n);
            if nargout >= 5
                Result_r = cell(1, n);
                Result_N = cell(1, n);
            end
            for t = n : -1 : 1
                if nargout >= 3
                    Result_r{t} = r;
                    Result_N{t} = N;
                end
                if Hdyn, H = Hmat{t}; end
                if Zdyn, Z = Zmat{t}; end
                if Tdyn, T = Tmat{t}; end
                if Qdyn, Q = Qmat{t}; end
                etahat(:, t) = QRt{t}*r;
                etavarhat(:,:, t) = Q - QRt{t}*N*RQ{t};
                if allmis(t)
                    %% Disturbance smoothing when all observations are missing %%
                    epshat(:, t) = 0;
                    epsvarhat(:,:, t) = H; %%%% TODO: What is epsvarhat when all missing?
                    r = T'*r;
                    N = T'*N*T;
                else
                    if anymis(t), Z(mis(:, t),:)=[]; H(mis(:, t),:)=[]; H(:,mis(:, t))=[]; end
                    if t > d || ~Fns(t)
                        %% Normal disturbance smoothing or when F_inf is zero %%
                        epshat(~mis(:, t), t) = H*(invF{t}*v{t} - K{t}'*r);
                        epsvarhat(~mis(:, t), ~mis(:, t), t) = H - H*(invF{t} + K{t}'*N*K{t})*H;
                        M = Z'*invF{t};
                        r = M*v{t} + L{t}'*r;
                        N = M*Z + L{t}'*N*L{t};
                    else
                        %% Exact initial disturbance smoothing when F_inf is nonsingular %%
                        epshat(~mis(:, t), t) = -H*K{t}'*r;
                        epsvarhat(~mis(:, t), ~mis(:, t), t) = H - H*K{t}'*N*K{t}*H;
                        r = L{t}'*r;
                        N = L{t}'*N*L{t};
                    end
                    if anymis(t), if ~Zdyn, Z = Zmat; end, if ~Hdyn, H = Hmat; end, end
                end
            end
        else [epshat etahat] = fastdisturbsmo_int(n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, ~issta(model.R), Qdyn, ~issta(model.c), Hmat, Zmat, Tmat, getmat(model.R), Qmat, getmat(model.c), a1, model.P1.mat, opt.tol);
        end
    end
    if nargout > 2
        %% Result postprocessing %%
        if size(epsvarhat, 1) == 1, epsvarhat = permute(epsvarhat, [2 3 1]); end
        if size(etavarhat, 1) == 1, etavarhat = permute(etavarhat, [2 3 1]); end
        if nargout > 4
            if ndims(y) > 2, Result_r = permute(cat(3, Result_r{:}), [1 3 2]); else Result_r = [Result_r{:}]; end
            Result_N    = cat(3, Result_N{:});
            if size(Result_N, 1) == 1, Result_N = permute(Result_N, [2 3 1]); end
        end
    end
end

