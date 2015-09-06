function [y alpha eps eta] = sample_int(n, N, p, m, r, Znl, Tnl, Z, T, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Rmat, Qmat, cmat, a1, P1)
% y is p*N*n
% alpha is m*N*n
% eps is p*N*n
% eta is r*N*n

%% Determine nonlinear functions %%
if Znl, Zdyn = false; else Zmat = getmat(Z); end
if Tnl, Tdyn = false; else Tmat = getmat(T); end

%% Draw from Gaussian distribution %%
alpha1  = randn(m, N);
eps     = randn(p, N, n);
eta     = randn(r, N, n);

%% Initialization for sampling %%
if Hdyn, eps(:,:, 1) = sigma_int(Hmat{1}, eps(:,:, 1));
else eps = reshape(sigma_int(Hmat, reshape(eps, p, N*n)), p, N, n); end
if Tdyn, T = Tmat{1}; elseif ~Tnl, T = Tmat; end
if Rdyn, R = Rmat{1}; else R = Rmat; end
if Qdyn, eta(:,:, 1) = sigma_int(Qmat{1}, eta(:,:, 1));
else eta = reshape(sigma_int(Qmat, reshape(eta, r, N*n)), r, N, n); end
if cdyn, c = repmat(cmat{1}, 1, N); else c = repmat(cmat, 1, N); end
P1(P1 == Inf) = 0;
y       = zeros(p, N, n);
alpha   = zeros(m, N, n);

%% Generate independent samples from the model %%
alpha(:,:, 1) = repmat(a1, 1, N) + sigma_int(P1, alpha1);
if Znl, y(:,:, 1) = getfunc(Z, alpha(:,:, 1), 1);
elseif Zdyn, y(:,:, 1) = Zmat{1}*alpha(:,:, 1);
end
for t = 2 : n
    if Hdyn, eps(:,:, t) = sigma_int(Hmat{t}, eps(:,:, t)); end
    if Qdyn, eta(:,:, t) = sigma_int(Qmat{t}, eta(:,:, t)); end
    if Tnl, alpha(:,:, t) = c + getfunc(T, alpha(:,:, t-1), t-1) + R*eta(:,:, t-1);
    else alpha(:,:, t) = c + T*alpha(:,:, t-1) + R*eta(:,:, t-1);
    end
    if Znl, y(:,:, t) = getfunc(Z, alpha(:,:, t), t); end
    if Zdyn, y(:,:, t) = Zmat{t}*alpha(:,:, t); end
    if Tdyn, T = Tmat{t}; end
    if Rdyn, R = Rmat{t}; end
    if cdyn, c = repmat(cmat{t}, 1, N); end
end
if ~Znl && ~Zdyn, y = reshape(Zmat*reshape(alpha, m, N*n), p, N, n) + eps;
else y = y + eps; end

%% Function for incorporating covariance into independent Gaussian samples %%
function x = sigma_int(Sigma, u)
dgSigma = diag(Sigma);
if isequal(Sigma, diag(dgSigma)), x = diag(sqrt(dgSigma))*u;
else % Sigma is not diagonal
    [U Lambda] = eig(full(Sigma));
    x = U*(diag(sqrt(diag(Lambda)))*u);
end
