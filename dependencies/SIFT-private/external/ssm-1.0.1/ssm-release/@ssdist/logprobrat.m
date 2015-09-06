function logw = logprobrat(d, N, varargin)

%@SSDIST/LOGPROBRAT Calculate the log probability ratio.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

%% Get simulated input %%
% inputs are all in x*n*N format
if nargin < 4
    eps     = varargin{1};
else
    y       = varargin{1};
    theta   = varargin{2};
    eps     = varargin{3};
    eps(isnan(eps)) = 0;
end

%% Calculate non-Gaussian log probabilities %%
logp    = zeros(1, N);
for i = 1 : N
    for j = 1 : length(d.type)
        diagmask    = d.diagmask(:, j);
        if d.type(j) % Additive noise
            logp(i) = logp(i) + sum(d.logpf{j}(eps(diagmask, :, i)));
        else % Exponential family
            logpj   = d.logpf{j}(y(diagmask, :), theta(diagmask, :, i));
            logpj(isnan(logpj)) = 0;
            logp(i) = logp(i) + sum(logpj);
        end
    end
end

%% Calculate Gaussian log probabilities %%
diagmask    = any(d.diagmask, 2);
if sum(d.diagmask, 1) == 1 % all non-Gaussian parts are restricted to diagonal elements
    ngmmask     = diag(diagmask);

    %% Preliminary calculations for loggeps %%
    H           = getdvec(d.ssmat, ngmmask);
    invH2       = 1./(2*H);
    logdetH2pi  = (reallog(H) + reallog(2*pi))/2;

    loggeps     = -repmat(invH2, [1 1 N]).*eps(diagmask, :, :).^2 - repmat(logdetH2pi, [1 1 N]);
    logg        = squeeze(sum(sum(loggeps, 1), 2))';
else
    %%%% TODO: Assumed that the Gaussian parts do not correlate with the non-Gaussian parts

    %% Preliminary calculations for loggeps %%
    n           = size(eps, 2);
    logdetH     = zeros(1, n);
    invH2       = cell(1, n);
    Hmat        = getmat(d.ssmat);
    for t = 1 : n
        H           = Hmat{t}(diagmask, diagmask);
        logdetH(t)  = reallog(det(H));
        invH2{t}    = inv(H)/2;
    end
    logdetH2pi  = (logdetH + size(H, 1)*reallog(2*pi))/2;

    loggeps     = zeros(N, n);
    eps         = permute(eps, [1 3 2]);
    for t = 1 : n
        loggeps(:, t)   = -diag(eps(:,:, t)'*invH2{t}*eps(:,:, t));
    end
    loggeps     = loggeps - repmat(logdetH2pi, N, 1);
    logg        = sum(loggeps, 2)';
end

%% Calculate log probability ratio %%
logw    = logp - logg;
