function [model logL output] = estimate(y, model, param0, alpha0, varargin)

%@SSMODEL/ESTIMATE Estimate the model parameters via MLE.
%   [model logL output] = ESTIMATE(y, model[, param0, alpha0])
%   ESTIMATE(..., optname1, optvalue1, optname2, optvalue2, ...)
%       y is the observation sequence.
%       model is the state space model to use.
%       param0 is the initial parameter values.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4 || isempty(alpha0), alpha0 = zeros(size(model.T, 1), size(y, 2));
else checkstate(model, alpha0, 'alpha0'); end
if nargin < 3 || isempty(param0), Xmask  = true(1, model.psi.w);
else
    if isnumeric(param0)
        if isscalar(param0), model.psi = set(model.psi, repmat(param0, 1, model.psi.w)); Xmask = true(1, model.psi.w);
        elseif size(param0, 2) ~= model.psi.w || ~any(size(param0, 1) == [1 2]), error('ssm:ssmodel:estimate:InputError', 'param0 must be a valid parameter value row vector for model.');
        elseif size(param0, 1) == 1, model.psi = set(model.psi, param0); Xmask = true(1, model.psi.w);
        else % size(param0, 1) == 2
            model.psi   = set(model.psi, param0(1, :));
            Xmask       = logical(param0(2, :));
        end
    elseif islogical(param0)
        if ~isequal(size(param0), [1 model.psi.w]), error('ssm:ssmodel:estimate:InputError', 'param0 must be a valid parameter value row vector for model.'); end
        Xmask   = param0;
    else error('ssm:ssmodel:estimate:InputError', 'param0 must be a valid parameter value row vector for model.');
    end
end

opt     = setopt([], varargin{:});
if strcmp(opt.fmin, 'bfgs') && ~exist('fminunc', 'file'), opt.fmin = 'simplex'; end

%% Data preprocessing %%
n       = size(y, 2);
mis     = isnan(y);
anymis  = any(mis, 1);
allmis  = all(mis, 1);
nmis    = n - nnz(allmis);

%% Get the current matrix values %%
psi     = model.psi.value;
p       = model.p;
m       = size(model.T, 1);
r       = size(model.R, 2);
w       = model.psi.w;
H       = model.H;
Z       = model.Z;
T       = model.T;
R       = model.R;
Q       = model.Q;
c       = model.c;
a1      = model.a1;
P1      = model.P1;
if opt.usec
    Hmat    = getmat_c(H);
    Zmat    = getmat_c(Z);
    Tmat    = getmat_c(T);
    Rmat    = getmat_c(R);
    Qmat    = getmat_c(Q);
    cmat    = getmat_c(c);
    a1mat   = getmat_c(a1);
    P1mat   = getmat_c(P1);
else
    Hmat    = getmat(H);
    Zmat    = getmat(Z);
    Tmat    = getmat(T);
    Rmat    = getmat(R);
    Qmat    = getmat(Q);
    cmat    = getmat(c);
    a1mat   = getmat(a1);
    P1mat   = getmat(P1);
end
Hdyn    = ~issta(H);
Zdyn    = ~issta(Z);
Tdyn    = ~issta(T);
Rdyn    = ~issta(R);
Qdyn    = ~issta(Q);
cdyn    = ~issta(c);
Hvar    = ~isconst(H);
Zvar    = ~isconst(Z);
Tvar    = ~isconst(T);
Rvar    = ~isconst(R);
Qvar    = ~isconst(Q);
cvar    = ~isconst(c);
a1var   = ~isconst(a1);
P1var   = ~isconst(P1);
Hdvar   = ~isdconst(H);
Zdvar   = ~isdconst(Z);
Tdvar   = ~isdconst(T);
Rdvar   = ~isdconst(R);
Qdvar   = ~isdconst(Q);
cdvar   = ~isdconst(c);
Hng     = isa(H, 'ssdist');
Qng     = isa(Q, 'ssdist');
Znl     = isa(Z, 'ssfunc');
Tnl     = isa(T, 'ssfunc');

%% Prepare data structure needed for model update %%
Ain             = model.A ~= 0;
Aout            = model.A > 0;
Azero           = ~any(Ain, 1); % We don't need these columns of A
Ain(:, Azero)   = [];
Aout(:, Azero)  = [];
adj_M           = logical(eye(18));
adj_M(:, Azero) = [];
adj_H           = find(adj_M(1, :));
adj_Z           = find(adj_M(2, :));
adj_T           = find(adj_M(3, :));
adj_R           = find(adj_M(4, :));
adj_Q           = find(adj_M(5, :));
adj_c           = find(adj_M(6, :));
adj_a1          = find(adj_M(7, :));
adj_P1          = find(adj_M(8, :));
adj_Hd          = find(adj_M(9, :));
adj_Zd          = find(adj_M(10, :));
adj_Td          = find(adj_M(11, :));
adj_Rd          = find(adj_M(12, :));
adj_Qd          = find(adj_M(13, :));
adj_cd          = find(adj_M(14, :));
adj_Hng         = find(adj_M(15, :));
adj_Qng         = find(adj_M(16, :));
adj_Znl         = find(adj_M(17, :));
adj_Tnl         = find(adj_M(18, :));
func            = model.func;
grad            = model.grad;
con             = model.psi.con;
pmask           = model.pmask;
if con, bounce  = confun(model.psi); end

%% Initialize search %%
if Hng || Qng
    % Non-Gaussian models
    % Get addition info needed for non-Gaussian case
    allHtype    = Hng && all(H.type);
    Hdistvar    = Hng && ~isdistconst(H);
    Qdistvar    = Qng && ~isdistconst(Q);
    N           = opt.nsamp;
    S           = randn('state');
    logw        = zeros(1, N);
    gradient    = false;
    fun         = @ngsnlogLfun;
elseif Znl || Tnl
    % Nonlinear models
    Zfvar       = Znl && ~isfconst(Z);
    Tfvar       = Tnl && ~isfconst(T);
    gradient    = false;
    fun         = @nlsnlogLfun;
else
    % Linear Gaussian models
    gradient    = ~Hdyn && ~Zvar && ~Zdvar && ~Tvar && ~Tdvar && ~Rvar && ~Rdyn && ~Qdyn ...
        && ~cvar && ~cdvar && ~a1var && ~P1var;
    if gradient
        for i = 1 : length(grad)
            if ~isa(grad{i}, 'function_handle'), gradient = false; break; end
        end
    end
    if gradient
        Hgrad   = zeros(size(Hmat));
        Qgrad   = zeros(size(Qmat));
        paramM  = logical(single(vertcat(pmask{:}))'*single(Aout));
        fun     = @snlogLfungrad;
    else
        fun     = @snlogLfun;
    end
end

%% Start searching %%
%%%% TODO: options for fmin?
switch opt.disp
    case 0, options = optimset('Display', 'off');
    case 1, options = optimset('Display', 'notify');
    case 2, options = optimset('Display', 'final'); opt.disp = 1; % No sub iterations shown
    case 3, options = optimset('Display', 'iter'); opt.disp = 1; % No sub iterations shown
end
switch opt.fmin
    case 'simplex'
        [psi1 L exitflag output]    = fminsearch(fun, psi(Xmask), options);
        output.exitflag             = exitflag;
    case 'bfgs'
        if gradient, options = optimset(options, 'GradObj', 'on', 'DerivativeCheck', 'off', 'LargeScale', 'off');
        else options = optimset(options, 'LargeScale', 'off');
        end
        [psi1 L exitflag output]    = fminunc(fun, psi(Xmask), options);
        output.exitflag             = exitflag;
    otherwise
        [psi1 L]    = feval(opt.fmin, fun, psi(Xmask));
end

%% Set model to maximum likelihood parameter values %%
psi(Xmask)  = psi1;
if con, psi = bounce(psi); end
model       = setparam(model, psi, true);
if Hng || Qng
    if opt.usec, [model.H model.Q output.ytilde] = gauss_int(n, y, [], [], [], Hng, Qng, model.H, model.Q, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, getmat_c(model.Z), getmat_c(model.T), getmat_c(model.R), getmat_c(model.c), getmat_c(model.a1), getmat_c(model.P1), alpha0, opt.tol, opt.maxiter, true, opt.inv);
    else [model.H model.Q output.ytilde] = gauss_int(n, y, mis, anymis, allmis, Hng, Qng, model.H, model.Q, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, getmat(model.Z), getmat(model.T), getmat(model.R), getmat(model.c), model.a1.mat, model.P1.mat, alpha0, opt.tol, opt.maxiter, false, opt.inv);
    end
elseif Znl || Tnl
    if opt.usec, [model.Z model.T model.c output.ytilde] = linear_int(n, y, [], [], [], Znl, Tnl, model.Z, model.T, model.c, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, getmat_c(model.H), getmat_c(model.R), getmat_c(model.Q), getmat_c(model.a1), getmat_c(model.P1), alpha0, opt.tol, opt.maxiter, true, opt.inv);
    else [model.Z model.T model.c output.ytilde] = linear_int(n, y, mis, anymis, allmis, Znl, Tnl, model.Z, model.T, model.c, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, getmat(model.H), getmat(model.R), getmat(model.Q), model.a1.mat, model.P1.mat, alpha0, opt.tol, opt.maxiter, false, opt.inv);
    end
end
logL        = -nmis*(p*reallog(2*pi) + L) / 2;
output.AIC  = (-2*logL + 2*(w + nnz(model.P1.mat == Inf)))/nmis;
output.BIC  = (-2*logL + reallog(nmis)*(w + nnz(model.P1.mat == Inf)))/nmis;

    %% Linear Gaussian function %%
    function L = snlogLfun(X)
        psi(Xmask) = X; X = psi;
        if con, X = bounce(X); end
        %% Update model %%
        Result  = cell(size(Ain));
        for i = 1 : length(func), [Result{i, Ain(i, :)}] = func{i}(X(pmask{i})); end
        if Hvar, H = setmat(H, vertcat(Result{Aout(:, adj_H), adj_H})); end
        if Zvar, Z = setmat(Z, vertcat(Result{Aout(:, adj_Z), adj_Z})); end
        if Tvar, T = setmat(T, vertcat(Result{Aout(:, adj_T), adj_T})); end
        if Rvar, R = setmat(R, vertcat(Result{Aout(:, adj_R), adj_R})); end
        if Qvar, Q = setmat(Q, vertcat(Result{Aout(:, adj_Q), adj_Q})); end
        if cvar, c = setmat(c, vertcat(Result{Aout(:, adj_c), adj_c})); end
        if Hdvar, H = setdvec(H, vertcat(Result{Aout(:, adj_Hd), adj_Hd})); end
        if Zdvar, Z = setdvec(Z, vertcat(Result{Aout(:, adj_Zd), adj_Zd})); end
        if Tdvar, T = setdvec(T, vertcat(Result{Aout(:, adj_Td), adj_Td})); end
        if Rdvar, R = setdvec(R, vertcat(Result{Aout(:, adj_Rd), adj_Rd})); end
        if Qdvar, Q = setdvec(Q, vertcat(Result{Aout(:, adj_Qd), adj_Qd})); end
        if cdvar, c = setdvec(c, vertcat(Result{Aout(:, adj_cd), adj_cd})); end
        if opt.usec
            if Hvar || Hdvar, Hmat = getmat_c(H); end
            if Zvar || Zdvar, Zmat = getmat_c(Z); end
            if Tvar || Tdvar, Tmat = getmat_c(T); end
            if Rvar || Rdvar, Rmat = getmat_c(R); end
            if Qvar || Qdvar, Qmat = getmat_c(Q); end
            if cvar || cdvar, cmat = getmat_c(c); end
            if a1var, a1mat = getmat_c(setmat(a1, vertcat(Result{Aout(:, adj_a1), adj_a1}))); end
            if P1var, P1mat = getmat_c(setmat(P1, vertcat(Result{Aout(:, adj_P1), adj_P1}))); end
            %% Calculate semi-negative loglikelihood %%
            L = loglik_int_c(y, Hmat, Hdyn, Zmat, Zdyn, Tmat, Tdyn, Rmat, Rdyn, Qmat, Qdyn, cmat, cdyn, a1mat, P1mat, opt.tol, opt.tol, opt.inv)/nmis;
        else
            if Hvar || Hdvar, Hmat = getmat(H); end
            if Zvar || Zdvar, Zmat = getmat(Z); end
            if Tvar || Tdvar, Tmat = getmat(T); end
            if Rvar || Rdvar, Rmat = getmat(R); end
            if Qvar || Qdvar, Qmat = getmat(Q); end
            if cvar || cdvar, cmat = getmat(c); end
            if a1var, a1mat = getmat(setmat(a1, vertcat(Result{Aout(:, adj_a1), adj_a1}))); end
            if P1var, P1mat = getmat(setmat(P1, vertcat(Result{Aout(:, adj_P1), adj_P1}))); end
            %% Calculate semi-negative loglikelihood %%
            L = kalman_int(4, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1mat, P1mat, opt.tol)/nmis;
        end
    end

    %% Linear Gaussian function and gradient %%
    function [f g] = snlogLfungrad(X)
        psi(Xmask) = X; X = psi;
        if con, X = bounce(X); end
        Result  = cell(size(Ain));
        if nargout < 2
            %% Update model %%
            for i = 1 : length(func), [Result{i, Ain(i, :)}] = func{i}(X(pmask{i})); end
            if opt.usec
                if Hvar, Hmat = getmat_c(setmat(H, vertcat(Result{Aout(:, adj_H), adj_H}))); end
                if Qvar, Qmat = getmat_c(setmat(Q, vertcat(Result{Aout(:, adj_Q), adj_Q}))); end
                %% Calculate semi-negative loglikelihood %%
                f = loglik_int_c(y, Hmat, Hdyn, Zmat, Zdyn, Tmat, Tdyn, Rmat, Rdyn, Qmat, Qdyn, cmat, cdyn, a1mat, P1mat, opt.tol, opt.tol, opt.inv)/nmis;
            else
                if Hvar, Hmat = getmat(setmat(H, vertcat(Result{Aout(:, adj_H), adj_H}))); end
                if Qvar, Qmat = getmat(setmat(Q, vertcat(Result{Aout(:, adj_Q), adj_Q}))); end
                %% Calculate semi-negative loglikelihood %%
                f = kalman_int(4, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1mat, P1mat, opt.tol)/nmis;
            end
        else
            %% Update model %%
            Result2 = cell(size(Ain));
            for i = 1 : length(grad), [Result{i, Ain(i, :)} Result2{i, Ain(i, :)}] = grad{i}(X(pmask{i})); end
            if opt.usec
                if Hvar, Hmat = getmat_c(setmat(H, vertcat(Result{Aout(:, adj_H), adj_H}))); end
                if Qvar, Qmat = getmat_c(setmat(Q, vertcat(Result{Aout(:, adj_Q), adj_Q}))); end
                [f rrN uuD] = loglik_grad_int_c(y, Hmat, Hdyn, Zmat, Zdyn, Tmat, Tdyn, Rmat, Rdyn, Qmat, Qdyn, cmat, cdyn, a1mat, P1mat, opt.tol, opt.tol, opt.inv);
                f = f/nmis;
                rrN = sum(rrN, 3); uuD = sum(uuD, 3);
                %% Calculate loglikelihood gradient %%
                Hgradvec    = blkdiag(Result2{Aout(:, adj_H), adj_H});
                Qgradvec    = blkdiag(Result2{Aout(:, adj_Q), adj_Q});
                g           = zeros(1, w);
                for i = 1 : w
                    if paramM(i, adj_H), Hgrad = getmat_c(setmat(H, Hgradvec(:, 1))); Hgradvec(:, 1) = []; else Hgrad(:, :) = 0; end
                    if paramM(i, adj_Q), Qgrad = getmat_c(setmat(Q, Qgradvec(:, 1))); Qgradvec(:, 1) = []; else Qgrad(:, :) = 0; end
                    g(i)    = (sum(diag(uuD*Hgrad)) + sum(diag(rrN*Rmat*Qgrad*Rmat')))/2;
                end
                g           = -2*g/nmis;
            else
                if Hvar, Hmat = getmat(setmat(H, vertcat(Result{Aout(:, adj_H), adj_H}))); end
                if Qvar, Qmat = getmat(setmat(Q, vertcat(Result{Aout(:, adj_Q), adj_Q}))); end
                %% Calculate semi-negative loglikelihood %%
                [d Fns v invF K L f] = kalman_int(8, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1mat, P1mat, opt.tol);
                f = f/nmis;
                %% Calculate smoothing cumulant %%
                [rrN uuD] = smocum_int(p, m, n, mis, anymis, allmis, Zdyn, Tdyn, Zmat, Tmat, d, Fns, v, invF, K, L);
                rrN = sum(rrN, 3); uuD = sum(uuD, 3);
                %% Calculate loglikelihood gradient %%
                Hgradvec    = blkdiag(Result2{Aout(:, adj_H), adj_H});
                Qgradvec    = blkdiag(Result2{Aout(:, adj_Q), adj_Q});
                g           = zeros(1, w);
                for i = 1 : w
                    if paramM(i, adj_H), Hgrad = getmat(setmat(H, Hgradvec(:, 1))); Hgradvec(:, 1) = []; else Hgrad(:, :) = 0; end
                    if paramM(i, adj_Q), Qgrad = getmat(setmat(Q, Qgradvec(:, 1))); Qgradvec(:, 1) = []; else Qgrad(:, :) = 0; end
                    g(i)    = (sum(diag(uuD*Hgrad)) + sum(diag(rrN*Rmat*Qgrad*Rmat')))/2;
                end
                g           = -2*g/nmis;
            end
        end
    end

    %% Non-Gaussian function %%
    function L = ngsnlogLfun(X)
        psi(Xmask) = X; X = psi;
        if con, X = bounce(X); end
        %% Update model %%
        Result  = cell(size(Ain));
        for i = 1 : length(func), [Result{i, Ain(i, :)}] = func{i}(X(pmask{i})); end
        if Hvar, H = setmat(H, vertcat(Result{Aout(:, adj_H), adj_H})); end
        if Zvar, Z = setmat(Z, vertcat(Result{Aout(:, adj_Z), adj_Z})); end
        if Tvar, T = setmat(T, vertcat(Result{Aout(:, adj_T), adj_T})); end
        if Rvar, R = setmat(R, vertcat(Result{Aout(:, adj_R), adj_R})); end
        if Qvar, Q = setmat(Q, vertcat(Result{Aout(:, adj_Q), adj_Q})); end
        if cvar, c = setmat(c, vertcat(Result{Aout(:, adj_c), adj_c})); end
        if Hdvar, H = setdvec(H, vertcat(Result{Aout(:, adj_Hd), adj_Hd})); end
        if Zdvar, Z = setdvec(Z, vertcat(Result{Aout(:, adj_Zd), adj_Zd})); end
        if Tdvar, T = setdvec(T, vertcat(Result{Aout(:, adj_Td), adj_Td})); end
        if Rdvar, R = setdvec(R, vertcat(Result{Aout(:, adj_Rd), adj_Rd})); end
        if Qdvar, Q = setdvec(Q, vertcat(Result{Aout(:, adj_Qd), adj_Qd})); end
        if cdvar, c = setdvec(c, vertcat(Result{Aout(:, adj_cd), adj_cd})); end
        if Hdistvar, H = setdist(H, vertcat(Result{Aout(:, adj_Hng), adj_Hng})); end
        if Qdistvar, Q = setdist(Q, vertcat(Result{Aout(:, adj_Qng), adj_Qng})); end

        if opt.usec
            if Zvar || Zdvar, Zmat = getmat_c(Z); end
            if Tvar || Tdvar, Tmat = getmat_c(T); end
            if Rvar || Rdvar, Rmat = getmat_c(R); end
            if cvar || cdvar, cmat = getmat_c(c); end
            if a1var, a1mat = getmat_c(setmat(a1, vertcat(Result{Aout(:, adj_a1), adj_a1}))); end
            if P1var, P1mat = getmat_c(setmat(P1, vertcat(Result{Aout(:, adj_P1), adj_P1}))); end
            %% Gaussian approximation %%
            [H Q ytilde]    = gauss_int(n, y, [], [], [], Hng, Qng, H, Q, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, ...
                Zmat, Tmat, Rmat, cmat, a1mat, P1mat, alpha0, opt.tol, opt.maxiter, true, opt.inv);
            Hmat            = getmat_c(H);
            Qmat            = getmat_c(Q);

            %% Calculate semi-negative loglikelihood %%
            L = loglik_int_c(ytilde, Hmat, Hdyn, Zmat, Zdyn, Tmat, Tdyn, Rmat, Rdyn, Qmat, Qdyn, cmat, cdyn, a1mat, P1mat, opt.tol, opt.tol, opt.inv);

            simsmo_int_c('seed', S);
            [alpha eps eta] = simsmo_int_c(ytilde, N, Hmat, Hdyn, Zmat, Zdyn, Tmat, Tdyn, Rmat, Rdyn, Qmat, Qdyn, cmat, cdyn, a1mat, P1mat, opt.antithetic, opt.tol, opt.tol, opt.inv, false);

            eps     = permute(eps, [1 3 2]);
            eta     = permute(eta, [1 3 2]);

            if ~allHtype
                if Zdyn, for t = 1:n, theta(:, t, :) = Zmat(:,:, t)*alpha(:,:, t); end
                else theta = permute(reshape(Zmat*reshape(alpha, 1, N*n), 1, N, n), [1 3 2]); end
            end
        else
            if Zvar || Zdvar, Zmat = getmat(Z); end
            if Tvar || Tdvar, Tmat = getmat(T); end
            if Rvar || Rdvar, Rmat = getmat(R); end
            if cvar || cdvar, cmat = getmat(c); end
            if a1var, a1mat = getmat(setmat(a1, vertcat(Result{Aout(:, adj_a1), adj_a1}))); end
            if P1var, P1mat = getmat(setmat(P1, vertcat(Result{Aout(:, adj_P1), adj_P1}))); end
            %% Gaussian approximation %%
            [H Q ytilde]    = gauss_int(n, y, mis, anymis, allmis, Hng, Qng, H, Q, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Zmat, Tmat, Rmat, cmat, a1mat, P1mat, alpha0, opt.tol, opt.maxiter, false, opt.inv);
            Hmat            = getmat(H);
            Qmat            = getmat(Q);

            %% Calculate Gaussian semi-negative logliklihood %%
            L = kalman_int(4, n, ytilde, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1mat, P1mat, opt.tol);

            %% Unconditional sampling %%
            randn('state', S);
            N   = 2*floor(N/2);
            [yplus alphaplus epsplus etaplus] = sample_int(n, N/2, p, m, r, Znl, Tnl, Z, T, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Rmat, Qmat, cmat, a1mat, P1mat);
            %% Generate conditional samples and antithetics %%
            if allHtype % alpha is not needed
                [epshat etahat] = fastdisturbsmo_int(n, ytilde, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1mat, P1mat, opt.tol);
                [epsplushat etaplushat] = batchdisturbsmo_int(n, N/2, yplus, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1mat, P1mat, opt.tol);
            else
                [alphahat epshat etahat] = fastsmo_int(n, ytilde, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1mat, P1mat, opt.tol);
                [alphaplushat epsplushat etaplushat] = batchsmo_int(n, N/2, yplus, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1mat, P1mat, opt.tol);
                alpha = repmat(permute(alphahat, [1 3 2]), [1 N 1]) + cat(2, -alphaplushat + alphaplus, alphaplushat - alphaplus);
                if Zdyn, for t = 1:n, theta(:, t, :) = Zmat{t}*alpha(:,:, t); end
                else theta = permute(reshape(Zmat*reshape(alpha, 1, N*n), 1, N, n), [1 3 2]); end
            end
            eps = repmat(epshat, [1 1 N]) + permute(cat(2, -epsplushat + epsplus, epsplushat - epsplus), [1 3 2]);
            eta = repmat(etahat, [1 1 N]) + permute(cat(2, -etaplushat + etaplus, etaplushat - etaplus), [1 3 2]);
        end

        %% Calculate non-Gaussian loglikelihood correction term %%
        if Hng
            if allHtype, logw = logprobrat(H, N, eps);
            else logw = logprobrat(H, N, y, theta, eps);
            end
        else logw = 0;
        end
        if Qng, logw = logw + logprobrat(Q, N, eta); end
        %% Calculate non-Gaussian semi-negative loglikelihood %%
        L = (L - 2*mean(logw))/nmis;
    end

    %% Nonlinear function %%
    function L = nlsnlogLfun(X)
        psi(Xmask) = X; X = psi;
        if con, X = bounce(X); end
        %% Update model %%
        Result  = cell(size(Ain));
        for i = 1 : length(func), [Result{i, Ain(i, :)}] = func{i}(X(pmask{i})); end
        if Hvar, H = setmat(H, vertcat(Result{Aout(:, adj_H), adj_H})); end
        if Zvar, Z = setmat(Z, vertcat(Result{Aout(:, adj_Z), adj_Z})); end
        if Tvar, T = setmat(T, vertcat(Result{Aout(:, adj_T), adj_T})); end
        if Rvar, R = setmat(R, vertcat(Result{Aout(:, adj_R), adj_R})); end
        if Qvar, Q = setmat(Q, vertcat(Result{Aout(:, adj_Q), adj_Q})); end
        if cvar, c = setmat(c, vertcat(Result{Aout(:, adj_c), adj_c})); end
        if a1var, a1 = setmat(a1, vertcat(Result{Aout(:, adj_a1), adj_a1})); end
        if P1var, P1 = setmat(P1, vertcat(Result{Aout(:, adj_P1), adj_P1})); end
        if Hdvar, H = setdvec(H, vertcat(Result{Aout(:, adj_Hd), adj_Hd})); end
        if Zdvar, Z = setdvec(Z, vertcat(Result{Aout(:, adj_Zd), adj_Zd})); end
        if Tdvar, T = setdvec(T, vertcat(Result{Aout(:, adj_Td), adj_Td})); end
        if Rdvar, R = setdvec(R, vertcat(Result{Aout(:, adj_Rd), adj_Rd})); end
        if Qdvar, Q = setdvec(Q, vertcat(Result{Aout(:, adj_Qd), adj_Qd})); end
        if cdvar, c = setdvec(c, vertcat(Result{Aout(:, adj_cd), adj_cd})); end
        if Zfvar, Z = setfunc(Z, vertcat(Result{Aout(:, adj_Znl), adj_Znl})); end
        if Tfvar, T = setfunc(T, vertcat(Result{Aout(:, adj_Tnl), adj_Tnl})); end

        if opt.usec
            if Hvar || Hdvar, Hmat = getmat_c(H); end
            if Rvar || Rdvar, Rmat = getmat_c(R); end
            if Qvar || Qdvar, Qmat = getmat_c(Q); end
            if a1var, a1mat = getmat_c(setmat(a1, vertcat(Result{Aout(:, adj_a1), adj_a1}))); end
            if P1var, P1mat = getmat_c(setmat(P1, vertcat(Result{Aout(:, adj_P1), adj_P1}))); end
            %% Linear approximation %%
            [Z T c ytilde]  = linear_int(n, y, [], [], [], Znl, Tnl, Z, T, c, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, ...
                Hmat, Rmat, Qmat, a1mat, P1mat, alpha0, opt.tol, opt.maxiter, true, opt.inv);
            Zmat            = getmat_c(Z);
            Tmat            = getmat_c(T);
            cmat            = getmat_c(c);
            cdyn            = Tnl;

            %% Calculate semi-negative loglikelihood %%
            L = loglik_int_c(ytilde, Hmat, Hdyn, Zmat, Zdyn, Tmat, Tdyn, Rmat, Rdyn, Qmat, Qdyn, cmat, cdyn, a1mat, P1mat, opt.tol, opt.tol, opt.inv)/nmis;
        else
            if Hvar || Hdvar, Hmat = getmat(H); end
            if Rvar || Rdvar, Rmat = getmat(R); end
            if Qvar || Qdvar, Qmat = getmat(Q); end
            if a1var, a1mat = getmat(setmat(a1, vertcat(Result{Aout(:, adj_a1), adj_a1}))); end
            if P1var, P1mat = getmat(setmat(P1, vertcat(Result{Aout(:, adj_P1), adj_P1}))); end
            %% Linear approximation %%
            [Z T c ytilde]  = linear_int(n, y, mis, anymis, allmis, Znl, Tnl, Z, T, c, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Rmat, Qmat, a1mat, P1mat, alpha0, opt.tol, opt.maxiter, false, opt.inv);
            Zmat            = getmat(Z);
            Tmat            = getmat(T);
            cmat            = getmat(c);
            cdyn            = Tnl;

            %% Calculate semi-negative logliklihood %%
            L = kalman_int(4, n, ytilde, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1mat, P1mat, opt.tol)/nmis;
        end
    end
end

