function [rrN uuD] = smocum_int(p, m, n, mis, anymis, allmis, Zdyn, Tdyn, Zmat, Tmat, d, Fns, v, invF, K, L)

%% Initialization %%
if ~Zdyn, Z = Zmat; end
if ~Tdyn, T = Tmat; end

%% Backwards recursion %%
r   = zeros(m, 1);
N   = zeros(m, m);
rrN = zeros(m, m, n);
uuD = zeros(p, p, n);
for t = n : -1 : 1
    rrN(:,:, t) = r*r' - N;
    if allmis(t)
        if Tdyn, T = Tmat{t}; end
        u   = zeros(p, 1);
        D   = zeros(p, p);
        r   = T'*r;
        N   = T'*N*T;
    else
        if t > d || ~Fns(t)
            if Zdyn, Z = Zmat{t}; end
            if anymis(t)
                Z(mis(:, t),:) = [];
                u   = zeros(p, 1);
                D   = zeros(p, p);
                u(~mis(:, t)) = invF{t}*v{t} - K{t}'*r;
                D(~mis(:, t), ~mis(:, t)) = invF{t} + K{t}'*N*K{t};
            else
                u   = invF{t}*v{t} - K{t}'*r;
                D   = invF{t} + K{t}'*N*K{t};
            end
            M   = Z'*invF{t};
            r   = M*v{t} + L{t}'*r;
            N   = M*Z + L{t}'*N*L{t};
            if anymis(t) && ~Zdyn, Z = Zmat; end
        else
            u   = -K{t}'*r;
            D   = K{t}'*N*K{t};
            r   = L{t}'*r;
            N   = L{t}'*N*L{t};
        end
    end
    uuD(:,:, t) = u*u' - D;
end
