function [w, V, invV, logdetV] = bayes_logit_fit_iter(X, y)
%% [w, V, invV, logdetV] = bayes_logit_fit_iter(X, y)
%
% returns parameters of fitted logit model
%
% p(y = 1 | x, w) = 1 / (1 + exp(- w' * x))
%
% with weight prior
%
% p(w) = N(w | 0, D^-1 I),
%
% where D is the size of x.
%
% The arguments are:
% X - input matrix, inputs x as row vectors
% y - output vector, containing either 1 or -1
%
% The function returns the posterior
%
% p(w1 | X, y) = N(w1 | w, V),
%
% and additionally the inverse of V and ln|V| (just in case).
%
% Compared to bayes_logit_fit, this function does not use a hyperprior,
% iterates over the inputs separately rather than processing them all at
% once, and is therefore slower, but also computationally more stable as
% it avoids computing the inverse of possibly close-to-singluar matrices.

% equations from Bishop (2006) PRML Book + errata (!)

[N, D] = size(X);

max_iter = 100;

%% (more or less) uninformative prior
V = eye(D) / D;
invV = eye(D) * D;
logdetV = - D * log(D);
w = zeros(D, 1);


%% iterate over all x separately
for n = 1:N;
    xn = X(n,:)';

    % precompute values
    Vx = V * xn;
    VxVx = Vx * Vx';
    c = xn' * Vx;
    xx = xn * xn';
    t_w = invV * w + 0.5 * y(n) * xn;

    % start iteration at xi = 0, lam_xi = 1/8
    V_xi = V - VxVx / (4 + c);
    invV_xi = invV + xx / 4;
    logdetV_xi = logdetV - log(1 + c / 4);
    w = V_xi * t_w;
    
    L_last = 0.5 * (logdetV_xi + w' * invV_xi * w) - log(2);

    for i = 1:max_iter
        % update xi by EM algorithm
        xi = sqrt(xn' * (V_xi + w * w') * xn);
        lam_xi = lam(xi);

        % Sherman-Morrison formula and Matrix determinant lemma
        V_xi = V - (2 * lam_xi / (1 + 2 * lam_xi * c)) * VxVx;
        invV_xi = invV + 2 * lam_xi * xx;
        logdetV_xi = logdetV - log(1 + 2 * lam_xi * c);
        w = V_xi * t_w;

        L = 0.5 * (logdetV_xi + w' * invV_xi * w - xi) ...
            - log(1 + exp(- xi)) + lam_xi * xi^2;

        % variational bound must grow!
        if L_last > L
            % fprintf('Last bound %6.6f, current bound %6.6f\n', L_last, L);
            % disp('Ran into numerical inaccuracies; terminating.');
            break;
        end
        % stop if change in variation bound is < 0.001%
        if abs(L_last - L) < abs(0.00001 * L)
            break
        end
        L_last = L;
    end
    if i == max_iter
        warning('Bayes:maxIter', ...
            'Bayesian logistic regression reached maximum number of iterations.');
    end
    
    V = V_xi;
    invV = invV_xi;
    logdetV = logdetV_xi;
end


function out = lam(xi)
% returns 1 / (4 * xi) * tanh(xi / 2)
divby0_w = warning('query', 'MATLAB:divideByZero');
warning('off', 'MATLAB:divideByZero');
out = tanh(xi ./ 2) ./ (4 .* xi);
warning(divby0_w.state, 'MATLAB:divideByZero');
% fix values where xi = 0
out(isnan(out)) = 1/8;
