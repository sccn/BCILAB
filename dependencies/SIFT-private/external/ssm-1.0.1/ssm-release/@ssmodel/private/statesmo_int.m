function varargout = statesmo_int(mode, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol)
% mode:
%   0 - all output.
%   1 - state smoother.
%   2 - ARMA EM. %%%% TODO: Ignore diffuse initialization for now.

%% Kalman filter %%
[a P d Fns v invF L P_inf F2 L1] = kalman_int(2, n, y, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Zmat, Tmat, Rmat, Qmat, cmat, a1, P1, tol);

%% Preallocate Output Results %%
switch mode
    case 0 % all output
        Result_r        = cell(1, n);
        Result_N        = cell(1, n);
    case 1 % state smoother
        Result_r        = cell(1, n);
        Result_N        = cell(1, n);
    case 2 % ARMA EM
        Result_N        = cell(1, n);
end

%% Initialization %%
if ~Zdyn, Z = Zmat; end
if ~Tdyn, T = Tmat; end

%% State smoothing backwards recursion %%
m   = size(a1, 1);
r   = zeros(m, 1);
r1  = zeros(m, 1);
N   = zeros(m, m);
N1  = zeros(m, m);
N2  = zeros(m, m);
alphahat    = zeros(size(a1, 1), n);
V           = zeros([size(P1) n]);
for t = n : -1 : 1
    switch mode
        case 0 % all output
            Result_r{t} = r;
            Result_N{t} = N;
        case 1 % state smoother
            Result_r{t} = r;
            Result_N{t} = N;
        case 2 % ARMA EM
            Result_N{t} = N;
    end
    if Zdyn, Z = Zmat{t}; end
    if Tdyn, T = Tmat{t}; end
    if allmis(t)
        %% State smoothing when all observations are missing %%
        r   = T'*r;
        N   = T'*N*T;
        if t > d
            alphahat(:, t)  = P(:,:, t)*r;
            V(:,:, t)       = P(:,:, t)*N*P(:,:, t);
        else
            r1  = T'*r1;
            N1  = T'*N1*T;
            N2  = T'*N2*T;
            alphahat(:, t)  = P(:,:, t)*r + P_inf(:,:, t)*r1;
            P_infN1P        = P_inf(:,:, t)*N1*P(:,:, t);
            V(:,:, t)       = P(:,:, t)*N*P(:,:, t) + P_infN1P' + P_infN1P + P_inf(:,:, t)*N2*P_inf(:,:, t);
        end
    else
        if anymis(t), Z(mis(:, t), :) = []; end
        if t > d
            %% Normal state smoothing %%
            M   = Z'*invF{t};
            r   = M*v{t} + L{t}'*r;
            N   = M*Z + L{t}'*N*L{t};
            alphahat(:, t)  = P(:,:, t)*r;
            V(:,:, t)       = P(:,:, t)*N*P(:,:, t);
        else
            %% Exact initial state smoothing %%
            if Fns(t)
                r1  = Z'*invF{t}*v{t} + L{t}'*r1 + L1{t}'*r;
                r   = L{t}'*r;
                LN  = L{t}'*N1 + L1{t}'*N;
                N2  = Z'*F2{t}*Z + LN*L1{t} + (L{t}'*N2 + L1{t}'*N1)*L{t};
                N1  = Z'*invF{t}*Z + LN*L{t};
                N   = L{t}'*N*L{t};
            else % F_inf{t} = 0
                M   = Z'*invF{t};
                r   = M*v{t} + L{t}'*r;
                r1  = T'*r1;
                N   = M*Z + L{t}'*N*L{t};
                N1  = T'*N1*L{t};
                N2  = T'*N2*T;
            end
            alphahat(:, t)  = P(:,:, t)*r + P_inf(:,:, t)*r1;
            P_infN1P        = P_inf(:,:, t)*N1*P(:,:, t);
            V(:,:, t)       = P(:,:, t)*N*P(:,:, t) + P_infN1P' + P_infN1P + P_inf(:,:, t)*N2*P_inf(:,:, t);
        end
        % Restore Z in case of stationarity and missing observations
        if anymis(t) && ~Zdyn, Z = Zmat; end
    end
end
alphahat    = a(:, 1:n) + alphahat;
V           = P(:,:, 1:n) - V;

%% Output Results %%
switch mode
    case 0 % all output
        varargout = {L P alphahat V Result_r Result_N};
    case 1 % state smoother
        varargout = {alphahat V Result_r Result_N};
    case 2 % ARMA EM
        varargout = {L P alphahat V Result_N};
end
