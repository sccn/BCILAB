function out = bayes_logit_post(X, w, V, invV)
%% out = bayes_logit_post(X, w, V, invV)
%
% returns a vector containing
%
% p(y=1 | x, X, Y)
%
% for x = each row in the given X, for a bayesian logit model with
%
% p(y = 1 | x, w1) = 1 / (1 + exp(- w1' * x)),
%
% and w, V, invV, logdetV are the posterior parameters N(w1 | w, V).

max_iter = 100;
[N, Dx] = size(X);

%% precompute some constants
w_t = ones(N, 1) * (w' * invV) + 0.5 * X; % w_t = V^-1 w + x / 2 as rows
Vx = X * V;                               % W x as rows
VxxVwt = Vx .* (sum(w_t .* Vx, 2) * ones(1, Dx)); % V x x^T V^T w_t as rows
Vwt = w_t * V;                            % V w_t as rows
xVx = sum(Vx .* X, 2);                    % x^T V x as rows
xVx2 = xVx .^ 2;                          % x^T V x x^T V x as rows

%% start first iteration with xi = 0, lam_xi = 1/8
xi = zeros(N, 1);
lam_xi = lam(xi);
a_xi = 1 ./ (4 + xVx);
w_xi = Vwt - (a_xi * ones(1, Dx)) .* VxxVwt;
logdetV_xi = - log(1 + xVx ./ 4);
wVw_xi = sum(w_xi .* (w_xi * invV), 2) + sum(w_xi .* X, 2) .^ 2 ./ 4;
L_last = 0.5 * (sum(logdetV_xi) + sum(wVw_xi)) - N * log(2);


%% iterate to from xi's that maximise variational bound
for i = 1:max_iter
    % update xi by EM algorithm
    xi = sqrt(xVx - a_xi .* xVx2 + sum(w_xi .* X, 2) .^ 2);
    lam_xi = lam(xi);

    % Sherman Morrison formula and Matrix determinant lemma
    a_xi = 2 .* lam_xi ./ (1 + 2 .* lam_xi .* xVx);
    w_xi = Vwt - (a_xi * ones(1, Dx)) .* VxxVwt;
    logdetV_xi = - log(1 + 2 * lam_xi .* xVx);

    % variational bound, omitting constant terms
    wVw_xi = sum(w_xi .* (w_xi * invV), 2) ...
        + 2 * lam_xi .* (sum(w_xi .* X, 2) .^ 2);
    L = sum(0.5 * (logdetV_xi + wVw_xi - xi) ...
            - log(1 + exp(-xi)) + lam_xi .* xi .^ 2);

    % variational bound must grow!
    if L_last > L
        fprintf('Last bound %6.6f, current bound %6.6f\n', L_last, L);
        break;
        %error('Variational bound should not reduce');
    end
    % stop if change in variation bound is < 0.001%
    if abs(L_last - L) < abs(0.00001 * L)
        break
    end
    L_last = L;
end


%% posterior from optimal xi's
out = 1 ./ (1 + exp(-xi)) ./ sqrt(1 + 2 .* lam_xi .* xVx) ...
    .* exp(0.5 .* (-xi - w' * invV * w + wVw_xi) + lam_xi .* xi .^ 2);


function out = lam(xi)
% returns 1 / (4 * xi) * tanh(xi / 2)
divby0_w = warning('query', 'MATLAB:divideByZero');
warning('off', 'MATLAB:divideByZero');
out = tanh(xi ./ 2) ./ (4 .* xi);
warning(divby0_w.state, 'MATLAB:divideByZero');
% fix values where xi = 0
out(isnan(out)) = 1/8;
