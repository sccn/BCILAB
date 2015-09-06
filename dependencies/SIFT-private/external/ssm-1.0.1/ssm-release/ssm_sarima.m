function model = ssm_sarima(p, d, q, P, D, Q, s, mean)

%SSM_SARIMA Create SSMODEL object for SARIMA(p, d, q)(P, D, Q)s model.
%   model = SSM_SARIMA(p, d, q, P, D, Q, s[, mean])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 7, error('ssm:ssm_sarima:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(p) || ~isnumeric(p), error('ssm:ssm_sarima:InputError', 'p must be a scalar.'); end
if ~isscalar(d) || ~isnumeric(d), error('ssm:ssm_sarima:InputError', 'd must be a scalar.'); end
if ~isscalar(q) || ~isnumeric(q), error('ssm:ssm_sarima:InputError', 'q must be a scalar.'); end
if ~isscalar(P) || ~isnumeric(P), error('ssm:ssm_sarima:InputError', 'P must be a scalar.'); end
if ~isscalar(D) || ~isnumeric(D), error('ssm:ssm_sarima:InputError', 'D must be a scalar.'); end
if ~isscalar(Q) || ~isnumeric(Q), error('ssm:ssm_sarima:InputError', 'Q must be a scalar.'); end
if ~isscalar(s) || ~isnumeric(s), error('ssm:ssm_sarima:InputError', 's must be a scalar.'); end
if nargin < 8, mean = false; end
[Z T R P1 Tmask Rmask P1mask]   = mat_sarima(p, d, q, P, D, Q, s, mean);
[fun gra psi]                   = fun_sarma(p, q, P, Q, s);
model                           = [ssm_null ssmodel(struct('type', 'sarima', 'p', p, 'd', d, 'q', q, 'P', P, 'D', D, 'Q', Q, 's', s), 0, Z, ssmat(T, Tmask), ssmat(R, Rmask), ssmat(0, true), [repmat('T', 1, p+P>0) repmat('R', 1, q+Q>0) 'QP1'], fun, gra, psi, [], ssmat(P1, P1mask))];

