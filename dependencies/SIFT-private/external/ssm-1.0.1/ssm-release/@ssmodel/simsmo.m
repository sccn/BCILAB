function [alphatilde epstilde etatilde alphaplus] = simsmo(y, model, N, varargin)

%@SSMODEL/SIMSMO Simulation smoothing with the linear Gaussian model.
%   [alphatilde epstilde etatilde] = SIMSMO(y, model, N)
%       generates N random state and disturbance sequences using model
%       conditional on data.
%   SIMSMO(..., optname1, optvalue1, optname2, optvalue2, ...)
%       y is the observation sequence.
%       model is the state space model to use.
%       N is the number of samples to generate.
%       alphatilde is the sampled states. (m*n*N)
%       epstilde is the sampled observation disturbances. (p*n*N)
%       etatilde is the sampled state disturbances. (m*n*N)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, N = 1; end

opt     = setopt([], varargin{:});

checkdata(model, y, 'y');

%%%% TODO: add more antithetic variables

Znl     = isa(model.Z, 'ssfunc');
Tnl     = isa(model.T, 'ssfunc');
Hdyn    = ~issta(model.H);
Zdyn    = ~issta(model.Z);
Tdyn    = ~issta(model.T);
Rdyn    = ~issta(model.R);
Qdyn    = ~issta(model.Q);
cdyn    = ~issta(model.c);

if opt.usec && ~Znl && ~Tnl
    [alphatilde epstilde etatilde alphaplus] = simsmo_int_c(y, N, getmat_c(model.H), Hdyn, getmat_c(model.Z), Zdyn, ...
        getmat_c(model.T), Tdyn, getmat_c(model.R), Rdyn, getmat_c(model.Q), Qdyn, getmat_c(model.c), cdyn, ...
        getmat_c(model.a1), getmat_c(model.P1), opt.antithetic, opt.tol, opt.tol, opt.inv, true);
else
    %% Get the current matrix values %%
    p       = model.p;
    m       = size(model.T, 1);
    r       = size(model.R, 2);
    Hmat    = getmat(model.H);
    Zmat    = getmat(model.Z);
    Tmat    = getmat(model.T);
    Rmat    = getmat(model.R);
    Qmat    = getmat(model.Q);
    cmat    = getmat(model.c);
    a1      = model.a1.mat;
    P1      = model.P1.mat;

    %% Data preprocessing %%
    n   = size(y, 2);
    mis = isnan(y);
    if opt.antithetic >= 1, Nsamp = ceil(N/2); else Nsamp = N; end

    %% Unconditional sampling %%
    [yplus alphaplus epsplus etaplus] = sample_int(n, Nsamp, p, m, r, Znl, Tnl, model.Z, model.T, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Rmat, Qmat, cmat, a1, P1);

    %% Fast (and batch) state and disturbance smoothing %%
    [alphahat epshat etahat] = fastsmo_int(n, y, mis, any(mis, 1), all(mis, 1), Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, opt.tol);
    [alphaplushat epsplushat etaplushat] = batchsmo_int(n, Nsamp, yplus, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, opt.tol);

    %% Calculate sampled disturbances, states and observations %%
    if opt.antithetic >= 1
        epstilde    = repmat(epshat, [1 1 2*Nsamp]) + permute(cat(2, -epsplushat + epsplus, epsplushat - epsplus), [1 3 2]);
        etatilde    = repmat(etahat, [1 1 2*Nsamp]) + permute(cat(2, -etaplushat + etaplus, etaplushat - etaplus), [1 3 2]);
        alphatilde  = repmat(alphahat, [1 1 2*Nsamp]) + permute(cat(2, -alphaplushat + alphaplus, alphaplushat - alphaplus), [1 3 2]);
        if(mod(N, 2))
            epstilde(:,:, N+1)      = [];
            etatilde(:,:, N+1)      = [];
            alphatilde(:,:, N+1)    = [];
        end
    else % antithetic <= 0
        epstilde    = repmat(epshat, [1 1 N]) - permute(epsplushat - epsplus, [1 3 2]);
        etatilde    = repmat(etahat, [1 1 N]) - permute(etaplushat - etaplus, [1 3 2]);
        alphatilde  = repmat(alphahat, [1 1 N]) - permute(alphaplushat - alphaplus, [1 3 2]);
    end
    alphaplus = permute(alphaplus, [1 3 2]);
end

