function [alphahat epshat etahat Result_r] = batchsmo_int(n, N, y, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol)
% y is p*N*n
% alphahat is m*N*n
% epshat is p*N*n
% etahat is r*N*n

%% Kalman filter %%
[d Fns v invF K L L1 QRt] = batchkalman_int(2, n, N, y, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol);

%% Initialization %%
a1      = repmat(a1, 1, N);
D       = (P1 == Inf);
P1(D)   = 0;
P1_inf  = double(D);
if ~Hdyn, H = Hmat; end
if ~Zdyn, Z = Zmat; end
if ~Tdyn, T = Tmat; end
if ~Rdyn, R = Rmat; end
if ~cdyn, c = repmat(cmat, 1, N); end

%% Disturbance smoothing backwards recursion %%
if nargout >= 4, Result_r = cell(1, n); end
m   = size(a1, 1);
r   = zeros(m, N);
r1  = zeros(m, N);
epshat  = zeros(size(y, 1), N, n);
etahat  = zeros(size(QRt{1}, 1), N, n);
for t = n : -1 : 1
    if nargout >= 4, Result_r{t} = r; end
    if Hdyn, H = Hmat{t}; end
    if Zdyn, Z = Zmat{t}; end
    etahat(:,:, t) = QRt{t}*r;
    if t > d || ~Fns(t)
        %% Normal disturbance smoothing or when F_inf is zero %%
        epshat(:,:, t) = H*(invF{t}*v{t} - K{t}'*r);
        r = Z'*invF{t}*v{t} + L{t}'*r;
        if t <= d
            if Tdyn, T = Tmat{t}; end
            r1 = T'*r1;
        end
    else
        %% Exact initial disturbance smoothing when F_inf is nonsingular %%
        epshat(:,:, t) = -H*K{t}'*r;
        r1 = Z'*invF{t}*v{t} + L{t}'*r1 + L1{t}'*r;
        r = L{t}'*r;
    end
end

%% Fast state smoothing %%
alphahat            = zeros(m, N, n);
alphahat(:,:, 1)    = a1 + P1*r + P1_inf*r1;
for t = 1 : n-1
    if Tdyn, T = Tmat{t}; end
    if Rdyn, R = Rmat{t}; end
    if cdyn, c = repmat(cmat{t}, 1, N); end
    alphahat(:,:, t+1) = c + T*alphahat(:,:, t) + R*etahat(:,:, t);
end
