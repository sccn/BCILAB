function dist = dist_multinomial(h, k)

%DIST_MULTINOMIAL Create SSDIST object for multinomial distribution.
%   dist = DIST_MULTINOMIAL(h, k)
%       h is the number of cells.
%       k is the number of trials at each time point, or a scalar if the
%           number of trials is stationary.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, error('ssm:dist_multinomial:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(h) || ~isscalar(h), error('ssm:dist_multinomial:InputError', 'h must be a scalar.'); end
if ~isnumeric(k) || size(k, 1) ~= 1, error('ssm:dist_multinomial:InputError', 'k must be a scalar or row vector.'); end

s_k         = k;
s_k1        = k + 1;
s_logfactk  = gammaln(k+1);
if ~isscalar(k), s_k2 = repmat(k, h-1, 1); else s_k2 = k; end

    function [H y] = matf_multinomial(y, theta)
        n       = size(theta, 2);
        Zdot    = exp(theta);
        B       = Zdot./repmat(sum(Zdot, 1) + 1, h-1, 1);
        db      = s_k2.*B;
        yb      = y - db;
        H       = cell(1, n);
        for t = 1 : n
            H{t}    = inv(diag(db(:, t)) - db(:, t)*B(:, t)');
            y(:, t) = H{t}*yb(:, t);
            H{t}    = H{t}(:);
        end
        H   = horzcat(H{:});
        y   = theta + y;
    end

    function logp = logpf_multinomial(y, theta)
        logp    = diag(y'*theta)' - s_k.*reallog(sum(exp(theta), 1) + 1) + s_logfactk - sum(gammaln(y+1), 1) - gammaln(s_k1 - sum(y, 1));
    end

dist    = ssdist(0, @matf_multinomial, @logpf_multinomial, size(k, 2), h-1);
end

