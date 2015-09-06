function [epshat etahat] = batchdisturbsmo_int(n, N, y, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol)
% y is p*N*n
% epshat is p*N*n
% etahat is r*N*n

%% Kalman filter %%
[d Fns v invF K L QRt] = batchkalman_int(4, n, N, y, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol);

%% Initialization %%
if ~Hdyn, H = Hmat; end
if ~Zdyn, Z = Zmat; end

%% Disturbance smoothing backwards recursion %%
m   = size(a1, 1);
r   = zeros(m, N);
epshat  = zeros(size(y, 1), N, n);
etahat  = zeros(size(QRt{1}, 1), N, n);
for t = n : -1 : 1
    if Hdyn, H = Hmat{t}; end
    etahat(:,:, t) = QRt{t}*r;
    if t > d || ~Fns(t)
        %% Normal disturbance smoothing or when F_inf is zero %%
        if Zdyn, Z = Zmat{t}; end
        epshat(:,:, t) = H*(invF{t}*v{t} - K{t}'*r);
        r = Z'*invF{t}*v{t} + L{t}'*r;
    else
        %% Exact initial disturbance smoothing when F_inf is nonsingular %%
        epshat(:,:, t) = -H*K{t}'*r;
        r = L{t}'*r;
    end
end
