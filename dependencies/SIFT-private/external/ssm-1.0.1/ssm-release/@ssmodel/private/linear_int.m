function [Z T c ytilde alpha converged iter] = linear_int(n, y, mis, anymis, allmis, Znl, Tnl, Z, T, c, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, Rmat, Qmat, a1, P1, alpha, tol, maxiter, usec, inv_method)

%% Linear approximation loop %%
ytilde      = y;
converged   = false;
iter        = maxiter + 1;
if Tnl
    cmmask  = any(T.horzmask, 2) > 0;
    %%%%%%%% TODO: proper support for existing c
    c       = ssmat(c.mat, [], cmmask);
    cdyn    = true;
end

if usec
    for i = 1 : maxiter
        if Znl
            [Z ytilde]  = setlinear(Z, alpha);
            ytilde      = y - sum(cat(3, ytilde{:}), 3);
        end
        if Tnl
            [T ctilde]  = setlinear(T, alpha);
            c           = setdvec(c, vertcat(ctilde{:}), cmmask);
        end
        if converged, break; end
        prevalpha   = alpha;
        alpha       = fastsmo_int_c(ytilde, Hmat, Hdyn, getmat_c(Z), Zdyn, getmat_c(T), Tdyn, Rmat, Rdyn, Qmat, Qdyn, ...
            getmat_c(c), cdyn, a1, P1, tol, tol, inv_method, true);
        if abs(alpha - prevalpha) < tol, converged = true; iter = i; end
    end
else
    for i = 1 : maxiter
        if Znl
            [Z ytilde]  = setlinear(Z, alpha);
            ytilde      = y - sum(cat(3, ytilde{:}), 3);
        end
        if Tnl
            [T ctilde]  = setlinear(T, alpha);
            c           = setdvec(c, vertcat(ctilde{:}), cmmask);
        end
        if converged, break; end
        prevalpha   = alpha;
        alpha       = faststatesmo_int(n, ytilde, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Hmat, getmat(Z), getmat(T), Rmat, Qmat, getmat(c), a1, P1, tol);
        if abs(alpha - prevalpha) < tol, converged = true; iter = i; end
    end
end
