function [epshat etahat] = fastdisturbsmo_int(n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol)

%% Kalman filter %%
[d Fns v invF K L QRt] = kalman_int(7, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol);

%% Initialization %%
if ~Hdyn, H = Hmat; end
if ~Zdyn, Z = Zmat; end
if ~Tdyn, T = Tmat; end

%% Disturbance smoothing backwards recursion %%
m   = size(a1, 1);
r   = zeros(m, 1);
epshat  = zeros(size(y, 1), n);
etahat  = zeros(size(QRt{1}, 1), n);
for t = n : -1 : 1
    etahat(:, t) = QRt{t}*r;
    if allmis(t)
        %% Disturbance smoothing when all observations are missing %%
        if Tdyn, T = Tmat{t}; end
        epshat(:, t) = 0;
        r = T'*r;
    else
        if Hdyn, H = Hmat{t}; end
        if Zdyn, Z = Zmat{t}; end
        if anymis(t), Z(mis(:, t),:)=[]; H(mis(:, t),:)=[]; H(:,mis(:, t))=[]; end
        if t > d || ~Fns(t)
            %% Normal disturbance smoothing or when F_inf is zero %%
            epshat(~mis(:, t), t) = H*(invF{t}*v{t} - K{t}'*r);
            r = Z'*invF{t}*v{t} + L{t}'*r;
        else
            %% Exact initial disturbance smoothing when F_inf is nonsingular %%
            epshat(~mis(:, t), t) = -H*K{t}'*r;
            r = L{t}'*r;
        end
        if anymis(t), if ~Zdyn, Z = Zmat; end, if ~Hdyn, H = Hmat; end, end
    end
end
