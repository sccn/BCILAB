function [logL fvar] = loglik(y, model, ytilde, varargin)

%@SSMODEL/LOGLIK Loglikelihood of data given model.
%   [logL fvar] = LOGLIK(y, model) calculates loglikelihood for the linear
%       Gaussian model.
%   logL = LOGLIK(y, model, ytilde) calculates loglikelihood for the
%       non-Gaussian model.
%   LOGLIK(..., optname1, optvalue1, optname2, optvalue2, ...)
%       y is the observation sequence.
%       model is the state space model to use.
%       logL is the loglikelihood of y.
%       fvar is the prediction error variance.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3 || isempty(ytilde), ytilde = y;
else checkdata(model, ytilde, 'ytilde'); end

opt     = setopt([], varargin{:});

checkdata(model, y, 'y');

if opt.usec
    %%---------------%%
    %% Gaussian part %%
    %%---------------%%

    %% Get the current matrix values %%
    Hd  = ~issta(model.H);
    Zd  = ~issta(model.Z);
    Td  = ~issta(model.T);
    Rd  = ~issta(model.R);
    Qd  = ~issta(model.Q);
    cd  = ~issta(model.c);
    H   = getmat_c(model.H);
    Z   = getmat_c(model.Z);
    T   = getmat_c(model.T);
    R   = getmat_c(model.R);
    Q   = getmat_c(model.Q);
    c   = getmat_c(model.c);
    a1  = getmat_c(model.a1);
    P1  = getmat_c(model.P1);

    %% Data preprocessing %%
    n       = size(ytilde, 2);
    mis     = isnan(ytilde);
    allmis  = all(mis, 1);

    %% Kalman filter %%
    [snlogL fvar] = loglik_int_c(ytilde, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, opt.tol, opt.tol, opt.inv);

    %% Calculate Gaussian loglikelihood %%
    logL    = -((n-nnz(allmis))*model.p*reallog(2*pi) + snlogL) / 2;
    fvar    = fvar / (n*model.p - nnz(P1 == Inf));

    %%-------------------%%
    %% Non-Gaussian part %%
    %%-------------------%%
    if isa(model.H, 'ssdist') || isa(model.Q, 'ssdist')
        %% Obtain simulation samples %%
        N   = opt.nsamp;
        [alpha eps eta] = simsmo_int_c(ytilde, N, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, ...
            opt.antithetic, opt.tol, opt.tol, opt.inv, false);

        if isa(model.H, 'ssdist')
            if ~all(model.H.type)
                for t = 1:n, theta(:, t, :) = Z(:,:, t)*alpha(:,:, t); end
                logw    = logprobrat(model.H, N, y, theta, eps);
            else logw   = logprobrat(model.H, N, eps);
            end
        end

        if isa(model.Q, 'ssdist'), logw = logw + logprobrat(model.Q, N, eta); end

        logL    = logL + mean(logw);
    end
else
    %%---------------%%
    %% Gaussian part %%
    %%---------------%%

    %% Get the current matrix values %%
    Zdyn    = ~issta(model.Z);
    Zmat    = getmat(model.Z);
    P1      = model.P1.mat;
    if ~Zdyn, Z = Zmat; end

    %% Data preprocessing %%
    n       = size(ytilde, 2);
    mis     = isnan(ytilde);
    allmis  = all(mis, 1);

    %% Kalman filter %%
    [logL_ fvar_] = kalman_int(4, n, ytilde, mis, any(mis, 1), allmis, ~issta(model.H), Zdyn, ~issta(model.T), ~issta(model.R), ~issta(model.Q), ~issta(model.c), getmat(model.H), Zmat, getmat(model.T), getmat(model.R), getmat(model.Q), getmat(model.c), model.a1.mat, P1, opt.tol);

    %% Calculate Gaussian loglikelihood %%
    logL    = -((n-nnz(allmis))*model.p*reallog(2*pi) + logL_) / 2;
    fvar    = fvar_ / (n*model.p - nnz(P1 == Inf));

    %%-------------------%%
    %% Non-Gaussian part %%
    %%-------------------%%
    if isa(model.H, 'ssdist') || isa(model.Q, 'ssdist')
        %% Obtain simulation samples %%
        N   = opt.nsamp;
        [alpha eps eta] = simsmo(ytilde, model, N, varargin{:});

        if isa(model.H, 'ssdist')
            if ~all(model.H.type)
                alpha   = permute(alpha, [1 3 2]);
                for t = 1 : n
                    if Zdyn, Z = Zmat{t}; end
                    theta(:, t, :) = Z*alpha(:, :, t);
                end
                logw    = logprobrat(model.H, N, y, theta, eps);
            else logw   = logprobrat(model.H, N, eps);
            end
        end

        if isa(model.Q, 'ssdist')
            logw    = logw + logprobrat(model.Q, N, eta);
        end

        logL    = logL + mean(logw);
    end
end

