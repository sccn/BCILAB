function varargout = diffdegree(y, varargin)

%DIFFDEGREE Determine differencing degree of data.
%   [d mean] = DIFFDEGREE(y[, ub, tsig])
%   [d D mean] = DIFFDEGREE(y, s[, ub, tsig])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2 || ~isscalar(varargin{1})
    % Unit root detection for ARIMA
    if nargin < 3, tsig = 1.5; else tsig = varargin{2}; end
    if nargin < 2, ub = [0.97 0.88]; else ub = varargin{1}; end
    arma1       = ssmodel('arma', 2, 0, true);
    arma1       = estimate(y, arma1, 0.1);
    phisq       = sqrt(arma1.param(1)^2 + 4*arma1.param(2));
    absphiroots = abs((arma1.param(1) + [phisq -phisq])/2);
    fprintf(1, 'Initial abs roots: %g, %g\n', absphiroots);
    if any(absphiroots > ub(1))
        d   = 1;
        y   = diff(y, d);
    else d  = 0;
    end
    arma2       = ssmodel('arma', 1, 1, true);
    while true
        arma2   = estimate(y, arma2, 0.1);
        fprintf(1, 'root: %g\n', arma2.param(1));
        if abs(arma2.param(1) - arma2.param(2)) < 0.15, break;
        elseif abs(arma2.param(1)) > ub(2)
            y   = diff(y);
            d   = d + 1;
        else break;
        end
    end
    [a P]       = kalman(y, arma2);
    tmean       = a(end, end)/realsqrt(P(end, end, end));
    fprintf(1, 'Final Mean: %g (%g)\n', a(end, end), tmean);
    varargout{1}    = d;
    varargout{2}    = abs(tmean) > tsig;
else
    % Unit root detection for SARIMA
    if nargin < 4, tsig = 1.5; else tsig = varargin{3}; end
    if nargin < 3, ub = [0.97 0.88]; else ub = varargin{2}; end
    s           = varargin{1};
    sarima1     = ssmodel('sarima', 2, 0, 0, 1, 0, 0, s, true);
    sarima1     = estimate(y, sarima1, 0.1);
    phisq       = sqrt(sarima1.param(1)^2 + 4*sarima1.param(2));
    absphiroots = abs((sarima1.param(1) + [phisq -phisq])/2);
    absPhiroot  = abs(sarima1.param(3));
    fprintf(1, 'Initial abs regular roots: %g, %g\n', absphiroots);
    fprintf(1, 'Initial abs seasonal root: %g\n', absPhiroot);
    if any(absphiroots > ub(1))
        d   = 1;
        y   = diff(y, d);
    else d  = 0;
    end
    if absPhiroot > ub(1)
        D   = 1;
        y   = y(s+1:end) - y(1:end-s);
    else D  = 0;
    end
    sarima2     = ssmodel('sarima', 1, 0, 1, 1, 0, 1, s, true);
    while true
        sarima2     = estimate(y, sarima2, 0.1);
        fprintf(1, 'regular root: %g\n', sarima2.param(1));
        fprintf(1, 'seasonal root: %g\n', sarima2.param(2));
        if abs(sarima2.param(1) - sarima2.param(3)) < 0.15, d2 = false;
        else d2 = abs(sarima2.param(1)) > ub(2); end
        if abs(sarima2.param(2) - sarima2.param(4)) < 0.15, D2 = false;
        else D2 = abs(sarima2.param(2)) > ub(2); end
        if ~d2 && ~D2, break;
        elseif d2 && D2 && d == 0 && D == 0
            if any(absphiroots >= absPhiroot), d = 1; y = diff(y);
            else D = 1; y = y(s+1:end) - y(1:end-s); end
        else
            if d2, d = d + 1; y = diff(y); end
            if D2, D = D + 1; y = y(s+1:end) - y(1:end-s); end
        end
    end
    [a P]       = kalman(y, sarima2);
    tmean       = a(end, end)/realsqrt(P(end, end, end));
    fprintf(1, 'Final Mean: %g (%g)\n', a(end, end), tmean);
    varargout{1}    = d;
    varargout{2}    = D;
    varargout{3}    = abs(tmean) > tsig;
end

