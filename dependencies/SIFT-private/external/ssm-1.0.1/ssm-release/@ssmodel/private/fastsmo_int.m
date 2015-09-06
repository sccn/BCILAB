function [alphahat epshat etahat] = fastsmo_int(n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol)

%% Kalman filter %%
[d Fns v invF K L L1 QRt] = kalman_int(5, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol);

%% Initialization %%
D       = (P1 == Inf);
P1(D)   = 0;
P1_inf  = double(D);
if ~Hdyn, H = Hmat; end
if ~Zdyn, Z = Zmat; end
if ~Tdyn, T = Tmat; end
if ~Rdyn, R = Rmat; end
if ~cdyn, c = cmat; end

%% Disturbance smoothing backwards recursion %%
m   = size(a1, 1);
r   = zeros(m, 1);
r1  = zeros(m, 1);
epshat  = zeros(size(y, 1), n);
etahat  = zeros(size(QRt{1}, 1), n);
for t = n : -1 : 1
    if Tdyn, T = Tmat{t}; end
    etahat(:, t) = QRt{t}*r;
    if allmis(t)
        %% Disturbance smoothing when all observations are missing %%
        epshat(:, t) = 0;
        r = T'*r;
        if t <= d, r1 = T'*r1; end
    else
        if Hdyn, H = Hmat{t}; end
        if Zdyn, Z = Zmat{t}; end
        if anymis(t), Z(mis(:, t),:)=[]; H(mis(:, t),:)=[]; H(:,mis(:, t))=[]; end
        if t > d || ~Fns(t)
            %% Normal disturbance smoothing or when F_inf is zero %%
            epshat(~mis(:, t), t) = H*(invF{t}*v{t} - K{t}'*r);
            r = Z'*invF{t}*v{t} + L{t}'*r;
            if t <= d, r1 = T'*r1; end
        else
            %% Exact initial disturbance smoothing when F_inf is nonsingular %%
            epshat(~mis(:, t), t) = -H*K{t}'*r;
            r1 = Z'*invF{t}*v{t} + L{t}'*r1 + L1{t}'*r;
            r = L{t}'*r;
        end
        if anymis(t), if ~Zdyn, Z = Zmat; end, if ~Hdyn, H = Hmat; end, end
    end
end

%% Fast state smoothing %%
alphahat        = zeros(m, n);
alphahat(:, 1)  = a1 + P1*r + P1_inf*r1;
for t = 1 : n-1
    if Tdyn, T = Tmat{t}; end
    if Rdyn, R = Rmat{t}; end
    if cdyn, c = cmat{t}; end
    alphahat(:, t+1) = c + T*alphahat(:, t) + R*etahat(:, t);
end
