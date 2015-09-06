function model = ssm_arima(p, d, q, mean)

%SSM_ARIMA Create SSMODEL object for ARIMA(p, d, q) models.
%   model = SSM_ARIMA(p, d, q[, mean])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, error('ssm:ssm_arima:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(p) || ~isnumeric(p), error('ssm:ssm_arima:InputError', 'p must be a scalar.'); end
if ~isscalar(d) || ~isnumeric(d), error('ssm:ssm_arima:InputError', 'd must be a scalar.'); end
if ~isscalar(q) || ~isnumeric(q), error('ssm:ssm_arima:InputError', 'q must be a scalar.'); end
if nargin < 4, mean = false; end
[Z T R P1 Tmask Rmask P1mask]   = mat_arima(p, d, q, mean);
[fun gra psi]                   = fun_arma(p, q);
model                           = [ssm_null ssmodel(struct('type', 'arima', 'p', p, 'd', d, 'q', q), 0, Z, ssmat(T, Tmask), ssmat(R, Rmask), ssmat(0, true), [repmat('T', 1, p>0) repmat('R', 1, q>0) 'QP1'], fun, gra, psi, [], ssmat(P1, P1mask))];

