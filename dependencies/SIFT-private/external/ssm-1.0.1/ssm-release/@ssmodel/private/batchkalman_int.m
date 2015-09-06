function varargout = batchkalman_int(mode, n, N, y, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a, P, tol)
% mode:
%   0 - all output.
%   1 - Kalman filter.
%   2 - fast smoother.
%   3 - fast state smoother.
%   4 - fast disturbance smoother.
%   y is p*N*n and has no missing data.
% (Only a and v depends on the data y, the values of P, P_inf, d, ... etc. is
% fixed for given model matrices (parameters).)

%% Initialization %%
a       = repmat(a, 1, N);
D       = (P == Inf);
init    = any(any(D)); % use exact initialization if init = true.
if init, d = n+1; P(D) = 0; P_inf = double(D); else d = 0; end
stationary  = ~Hdyn && ~Zdyn && ~Tdyn && ~Rdyn && ~Qdyn;
converged   = false;
RQdyn       = Rdyn || Qdyn;
if ~Hdyn, H = Hmat; end
if ~Zdyn, Z = Zmat; end
if ~Tdyn, T = Tmat; end
if ~Rdyn, R = Rmat; end
if ~Qdyn, Q = Qmat; end
if ~cdyn, c = repmat(cmat, 1, N); end
if ~RQdyn, RQRt = R*Q*R'; end

%% Preallocate Output Results %%
switch mode
    case 0 % all output
        Result_a = zeros(size(a, 1), N, n+1);
        Result_a(:,:, 1) = a;
        Result_P = zeros([size(P) n+1]);
        Result_P(:,:, 1) = P;
        Result_v    = cell(1, n);
        Result_invF = cell(1, n);
        Result_K    = cell(1, n);
        Result_L    = cell(1, n);
        Result_Pinf = zeros([size(P) nnz(D)]);
        if init, Result_Pinf(:,:, 1) = P_inf; end
        Result_F2   = cell(1, nnz(D));
        Result_L1   = cell(1, nnz(D));
        Result_RQ   = cell(1, n);
        if ~RQdyn, [Result_RQ{:}] = deal(R*Q); end
        Result_QRt  = cell(1, n);
        if ~RQdyn, [Result_QRt{:}] = deal(Q*R'); end
        Result_RQRt = cell(1, n);
        if ~RQdyn, [Result_RQRt{:}] = deal(RQRt); end
    case 1 % Kalman filter
        Result_a = zeros(size(a, 1), N, n+1);
        Result_a(:,:, 1) = a;
        Result_P = zeros([size(P) n+1]);
        Result_P(:,:, 1) = P;
        Result_v    = cell(1, n);
        Result_invF = cell(1, n);
    case 2 % fast smoother
        Result_v    = cell(1, n);
        Result_invF = cell(1, n);
        Result_K    = cell(1, n);
        Result_L    = cell(1, n);
        Result_L1   = cell(1, nnz(D));
        Result_QRt  = cell(1, n);
        if ~RQdyn, [Result_QRt{:}] = deal(Q*R'); end
    case 3 % fast state smoother
        Result_v    = cell(1, n);
        Result_invF = cell(1, n);
        Result_L    = cell(1, n);
        Result_L1   = cell(1, nnz(D));
        Result_RQRt = cell(1, n);
        if ~RQdyn, [Result_RQRt{:}] = deal(RQRt); end
    case 4 % fast disturbance smoother
        Result_v    = cell(1, n);
        Result_invF = cell(1, n);
        Result_K    = cell(1, n);
        Result_L    = cell(1, n);
        Result_QRt  = cell(1, n);
        if ~RQdyn, [Result_QRt{:}] = deal(Q*R'); end
end
Fns = true(1, n); % Is F_inf nonsingular for each iteration

%% Batch Kalman filter loop %%
for t = 1 : n
    if Hdyn, H = Hmat{t}; end
    if Zdyn, Z = Zmat{t}; end
    if Tdyn, T = Tmat{t}; end
    if Rdyn, R = Rmat{t}; end
    if Qdyn, Q = Qmat{t}; end
    if cdyn, c = repmat(cmat{t}, 1, N); end
    if RQdyn, RQRt = R*Q*R'; end
    if ~converged
        if init
            %% Exact initial Kalman filter %%
            M       = P*Z';
            M_inf   = P_inf*Z';
            A_inf   = T*P_inf;
            if abs(M_inf) < tol % F_inf is zero
                Fns(t)  = false;
                F       = Z*M + H;
                F2(:)   = 0;
                invF    = inv(F);
                K       = T*M*invF;
                K1(:)   = 0;
                L       = T - K*Z;
                L1(:)   = 0;
                P       = T*P*L' + RQRt;
                P_inf   = A_inf*T';
            else % F_inf is assumed to be nonsingular
                invF    = inv(Z*M_inf);
                F       = Z*M + H;
                F2      = -invF*F*invF;
                K       = T*M_inf*invF;
                K1      = T*(M*invF + M_inf*F2);
                L       = T - K*Z;
                L1      = -K1*Z;
                P       = T*P*L' + A_inf*L1' + RQRt;
                P_inf   = A_inf*L';
            end
            if abs(P_inf) < tol, d = t; init = false; end
        else
            %% Normal Kalman filter %%
            M       = P*Z';
            F       = Z*M + H;
            invF    = inv(F);
            K       = T*M*invF;
            L       = T - K*Z;
            prevP   = P;
            P       = T*P*L' + RQRt;
            if stationary, if abs(P-prevP) < tol, converged = true; end, end
        end
    end
    %% Kalman data filter %%
    v   = y(:,:, t) - Z*a;
    a   = c + T*a + K*v;

    %% Store results for this time point %%
    switch mode
        case 0 % all output
            Result_a(:,:, t+1) = a;
            Result_P(:,:, t+1) = P;
            Result_v{t}    = v;
            Result_invF{t} = invF;
            Result_K{t}    = K;
            Result_L{t}    = L;
            if t <= d
                Result_Pinf(:,:, t+1) = P_inf;
                if Fns(t)
                    Result_F2{t} = F2;
                    Result_L1{t} = L1;
                end
            end
            if RQdyn
                Result_RQ{t}   = R*Q;
                Result_QRt{t}  = Q*R';
                Result_RQRt{t} = RQRt;
            end
        case 1 % Kalman filter
            Result_a(:,:, t+1) = a;
            Result_P(:,:, t+1) = P;
            Result_v{t}    = v;
            Result_invF{t} = invF;
        case 2 % fast smoother
            Result_v{t}    = v;
            Result_invF{t} = invF;
            Result_K{t}    = K;
            Result_L{t}    = L;
            if t <= d && Fns(t), Result_L1{t} = L1; end
            if RQdyn, Result_QRt{t}  = Q*R'; end
        case 3 % fast state smoother
            Result_v{t}    = v;
            Result_invF{t} = invF;
            Result_L{t}    = L;
            if t <= d && Fns(t), Result_L1{t} = L1; end
            if RQdyn, Result_RQRt{t} = RQRt; end
        case 4 % fast disturbance smoother
            Result_v{t}    = v;
            Result_invF{t} = invF;
            Result_K{t}    = K;
            Result_L{t}    = L;
            if RQdyn, Result_QRt{t}  = Q*R'; end
    end
end
%% Output Results %%
switch mode
    case 0 % all output
        varargout = {Result_a Result_P d Fns Result_v Result_invF Result_K Result_L Result_Pinf Result_F2 Result_L1 Result_RQ Result_QRt Result_RQRt};
    case 1 % Kalman filter
        varargout = {Result_a Result_P d Result_v Result_invF};
    case 2 % fast smoother
        varargout = {d Fns Result_v Result_invF Result_K Result_L Result_L1 Result_QRt};
    case 3 % fast state smoother
        varargout = {d Fns Result_v Result_invF Result_L Result_L1 Result_RQRt};
    case 4 % fast disturbance smoother
        varargout = {d Fns Result_v Result_invF Result_K Result_L Result_QRt};
end
