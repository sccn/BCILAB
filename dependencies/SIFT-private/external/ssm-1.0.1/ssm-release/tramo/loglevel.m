function lambda = loglevel(y, varargin)

%LOGLEVEL Test data for the log-level specification.
%   lambda = LOGLEVEL(y, s)
%   lambda = LOGLEVEL(y, p, q)
%   lambda = LOGLEVEL(y, p, d, q)
%   lambda = LOGLEVEL(y, p, d, q, P, D, Q, s)
%       y is the observation data.
%       s is the seasonal period.
%       lambda is 1 if no logs are needed, else 0.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

switch nargin
    case 2
        if isa(varargin{1}, 'ssmodel'), model = varargin{1};
        else model = ssmodel('airline', varargin{1}); end
    case 3
        model   = ssmodel('arma', varargin{:});
    case 4
        model   = ssmodel('arima', varargin{:});
    case 8
        model   = ssmodel('sarima', varargin{:});
    otherwise
        model   = ssmodel('airline');
end
[model logL0]   = estimate(reallog(y), model, 0.1);
logL0           = logL0 - sum(reallog(y)); % Correction for Jacobian of log transform
[model logL1]   = estimate(y, model, 0.1);
fprintf(1, 'logL0 = %g, logL1 = %g\n', logL0, logL1);
lambda          = logL0 <= logL1;

