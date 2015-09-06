function alphahat = faststatesmo_int(n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol)

%% Kalman filter %%
[d Fns v invF L L1 RQRt] = kalman_int(6, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol);

%% Initialization %%
D       = (P1 == Inf);
P1(D)   = 0;
P1_inf  = double(D);
if ~Zdyn, Z = Zmat; end
if ~Tdyn, T = Tmat; end
if ~cdyn, c = cmat; end

%% Backwards recursion %%
m   = size(a1, 1);
r   = zeros(m, n+1);
r1  = zeros(m, 1);
for t = n : -1 : 1
    if Tdyn, T = Tmat{t}; end
    if allmis(t)
        %% Disturbance smoothing when all observations are missing %%
        r(:, t) = T'*r(:, t+1);
        if t <= d, r1 = T'*r1; end
    else
        if Zdyn, Z = Zmat{t}; end
        if anymis(t), Z(mis(:, t),:)=[]; end
        if t > d || ~Fns(t)
            %% Normal disturbance smoothing or when F_inf is zero %%
            r(:, t) = Z'*invF{t}*v{t} + L{t}'*r(:, t+1);
            if t <= d, r1 = T'*r1; end
        else
            %% Exact initial disturbance smoothing when F_inf is nonsingular %%
            r(:, t) = L{t}'*r(:, t+1);
            r1 = Z'*invF{t}*v{t} + L{t}'*r1 + L1{t}'*r(:, t+1);
        end
        if anymis(t) && ~Zdyn, Z = Zmat; end
    end
end

%% Fast state smoothing %%
alphahat        = zeros(m, n);
alphahat(:, 1)  = a1 + P1*r(:, 1) + P1_inf*r1;
for t = 1 : n-1
    if Tdyn, T = Tmat{t}; end
    if cdyn, c = cmat{t}; end
    alphahat(:, t+1) = c + T*alphahat(:, t) + RQRt{t}*r(:, t+1);
end
