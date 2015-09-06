function varargout = weights_int(mode, n, t0, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, Hmat, Zmat, Tmat, Rmat, Qmat, P, tol)
% mode:
%   0 - all primitive outputs.
%   1 - calculate Kalman filter weights.
%   2 - calculate Kalman smoother weights.

%% Initialization %%
m       = size(P, 1);
D       = (P == Inf);
init    = any(any(D)); % use exact initialization if init = true.
if init, d = n+1; P(D) = 0; P_inf = double(D);
else d = 0; end
stationary  = ~Hdyn && ~Zdyn && ~Tdyn && ~Rdyn && ~Qdyn;
converged   = false;
RQdyn       = Rdyn || Qdyn;
if ~Hdyn, H = Hmat; end
if ~Zdyn, Z = Zmat; end
if ~Tdyn, T = Tmat; end
if ~Rdyn, R = Rmat; end
if ~Qdyn, Q = Qmat; end
if ~RQdyn, RQRt = R*Q*R'; end

%% Preallocate Output Results %%
switch mode
    case 0 % all primitive outputs
        Result_P    = cell(1, n);
        Result_P{1} = P;
        Result_invF = cell(1, n);
        Result_K    = cell(1, n);
        Result_L    = cell(1, n);
        Result_N    = cell(1, n);
        Result_W    = cell(1, n);
    case 1 % calculate Kalman filter weights
        Result_K        = cell(1, t0-1);
        Result_L        = cell(1, t0-1);
        Result_omega    = cell(1, t0-1);
    case 2 % calculate Kalman smoother weights
        Result_invF         = cell(1, n);
        Result_K            = cell(1, n);
        Result_L            = cell(1, n);
        Result_LLasc        = cell(1, n);
        if t0 > 1, Result_LLasc{t0-1} = eye(m);
        else Result_P = P; end % t0 == 1
        Result_omega        = cell(1, t0-1);
        Result_omegaalpha   = cell(1, n);
end
Fns     = true(1, n); % Is F_inf nonsingular for each iteration

%% Kalman filter loop %%
for t = 1 : n
    if Hdyn, H = Hmat{t}; end
    if Zdyn, Z = Zmat{t}; end
    if Tdyn, T = Tmat{t}; end
    if Rdyn, R = Rmat{t}; end
    if Qdyn, Q = Qmat{t}; end
    if RQdyn, RQRt = R*Q*R'; end
    if ~converged
        if init
            %% Exact initial Kalman filter %%
            M       = P*Z';
            M_inf   = P_inf*Z';
            A_inf   = T*P_inf;
            if abs(M_inf) < tol % F_inf is zero
                Fns(t)  = false;
                invF    = inv(Z*M + H);
                K       = T*M*invF;
                L       = T - K*Z;
                P       = T*P*L' + RQRt;
                P_inf   = A_inf*T';
            else % F_inf is assumed to be nonsingular
                invF    = inv(Z*M_inf); % This is actually invF1
                F2      = -invF*(Z*M + H)*invF;
                K       = T*M_inf*invF;
                K1      = T*(M*invF + M_inf*F2);
                L       = T - K*Z;
                L1      = -K1*Z;
                P       = A_inf*L1' + T*P*L' + RQRt;
                P_inf   = A_inf*L';
            end
            if abs(P_inf) < tol, d=t; init=false; end
        else
            %% Normal Kalman filter %%
            M       = P*Z';
            invF    = inv(Z*M + H);
            K       = T*M*invF;
            L       = T - K*Z;
            prevP   = P;
            P       = T*P*L' + RQRt;
            if stationary, if abs(P-prevP) < tol, converged = true; end, end
        end
    end
    %% Store results for this time point %%
    switch mode
        case 0 % all primitive outputs
            Result_P{t+1}   = P;
            Result_invF{t}  = invF;
            Result_K{t}     = K;
            Result_L{t}     = L;
        case 1 % calculate Kalman filter weights
            Result_K{t}     = K;
            Result_L{t}     = L;
            if t >= t0 - 1, break; end
        case 2 % calculate Kalman smoother weights
            if t == t0 - 1, Result_P    = P; end
            if t == t0, Result_LLasc{t} = L'; end
            Result_invF{t}  = invF;
            Result_K{t}     = K;
            Result_L{t}     = L;
            if t > t0
                Result_LLasc{t} = Result_LLasc{t-1}*L';
            end
    end
end

%% Backwards recursion %%
switch mode
    case 0 % all primitive outputs
        N       = zeros(m);
        for t = n : -1 : 1
            if Hdyn, H = Hmat{t}; end
            if Zdyn, Z = Zmat{t}; end
            invF        = Result_invF{t};
            K           = Result_K{t};
            L           = Result_L{t};
            Result_W{t} = H*(invF*Z - K'*N*L);
            Result_N{t} = N;
            if t <= d && Fns(t), N = L'*N*L; else N = Z'*invF*Z + L'*N*L; end
        end
    case 1 % calculate Kalman filter weights
        LLdsc   = eye(m);
        for t = t0-1 : -1 : 1
            Result_omega{t} = LLdsc*Result_K{t};
            LLdsc           = LLdsc*Result_L{t};
        end
    case 2 % calculate Kalman smoother weights
        N       = zeros(m);
        LLdsc   = eye(m);
        for t = n : -1 : 1
            if t >= t0
                if t <= d && Fns(t)
                    K       = Result_K{t};
                    L       = Result_L{t};
                    W2      = -K'*N*L;
                    N       = L'*N*L;
                else
                    if Zdyn, Z = Zmat{t}; end
                    invF    = Result_invF{t};
                    K       = Result_K{t};
                    L       = Result_L{t};
                    W2      = invF*Z - K'*N*L;
                    N       = Z'*invF*Z + L'*N*L;
                end
                if t == t0
                    Result_omegaalpha{t} = Result_P*W2';
                    IPN = eye(m) - Result_P*N;
                else Result_omegaalpha{t} = Result_P*Result_LLasc{t-1}*W2'; end
            else
                Result_omega{t}         = LLdsc*Result_K{t};
                LLdsc                   = LLdsc*Result_L{t};
                Result_omegaalpha{t}    = IPN*Result_omega{t};
            end
        end
end

%% Output Results %%
switch mode
    case 0 % all primitive outputs
        varargout = {Result_P Result_invF Result_K Result_L Result_W Result_N};
    case 1 % calculate Kalman filter weights
        varargout = {Result_omega};
    case 2 % calculate Kalman smoother weights
        varargout = {Result_omega Result_omegaalpha};
end
