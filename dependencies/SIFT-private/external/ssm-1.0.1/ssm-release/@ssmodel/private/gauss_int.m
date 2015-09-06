function [H Q ytilde alpha converged c] = gauss_int(n, y, mis, anymis, allmis, Hng, Qng, H, Q, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, Zmat, Tmat, Rmat, cmat, a1, P1, alpha, tol, maxiter, usec, inv_method)

%% Gaussian approximation loop %%
ytilde      = y;
converged   = false;
c           = maxiter + 1;

if usec
    for i = 1 : maxiter
        if Hng
            if Zdyn, for t = 1:n, theta(:, t) = Zmat(:,:, t)*alpha(:, t); end
            else theta = Zmat*alpha; end
            [H ytilde] = setgauss(H, y, theta);
        end
        if Qng
            if Tdyn, for t = 1:n-1, eta(:, t) = alpha(:, t+1) - Tmat(:,:, t)*alpha(:, t); end
            else eta = alpha(:, 2:end) - Tmat*alpha(:, 1:end-1); end
            eta(:, n) = 0;
            Q = setgauss(Q, eta);
        end
        if converged, break; end
        prevalpha   = alpha;
        alpha       = fastsmo_int_c(ytilde, getmat_c(H), Hdyn, Zmat, Zdyn, Tmat, Tdyn, Rmat, Rdyn, getmat_c(Q), Qdyn, ...
            cmat, cdyn, a1, P1, tol, tol, inv_method, true);
        if abs(alpha - prevalpha) < tol, converged = true; c = i; end
    end
else
    for i = 1 : maxiter
        if Hng
            if Zdyn, for t = 1:n, theta(:, t) = Zmat{t}*alpha(:, t); end
            else theta = Zmat*alpha; end
            [H ytilde] = setgauss(H, y, theta);
        end
        if Qng
            if Tdyn, for t = 1:n-1, eta(:, t) = alpha(:, t+1) - Tmat{t}*alpha(:, t); end
            else eta = alpha(:, 2:end) - Tmat*alpha(:, 1:end-1); end
            eta(:, n) = 0;
            Q = setgauss(Q, eta);
        end
        if converged, break; end
        prevalpha   = alpha;
        alpha = faststatesmo_int(n, ytilde, mis, anymis, allmis, Hdyn, Zdyn, Tdyn, Rdyn, Qdyn, cdyn, getmat(H), Zmat, Tmat, Rmat, getmat(Q), cmat, a1, P1, tol);
        if abs(alpha - prevalpha) < tol, converged = true; c = i; end
    end
end
