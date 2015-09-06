function model = ssm_sarimahtd(p, d, q, P, D, Q, s, gauss)

%SSM_SARIMAHTD Create SSMODEL object for SARIMA with Hillmer-Tiao Decomposition and non-Gaussian noise.
%   model = SSM_SARIMAHTD(p, d, q, P, D, Q, s[, gauss])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 7, error('ssm:ssm_sarimahtd:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(p) || ~isnumeric(p), error('ssm:ssm_sarimahtd:InputError', 'p must be a scalar.'); end
if ~isscalar(d) || ~isnumeric(d), error('ssm:ssm_sarimahtd:InputError', 'd must be a scalar.'); end
if ~isscalar(q) || ~isnumeric(q), error('ssm:ssm_sarimahtd:InputError', 'q must be a scalar.'); end
if ~isscalar(P) || ~isnumeric(P), error('ssm:ssm_sarimahtd:InputError', 'P must be a scalar.'); end
if ~isscalar(D) || ~isnumeric(D) || D < 1, error('ssm:ssm_sarimahtd:InputError', 'D must be a scalar no less than one.'); end
if ~isscalar(Q) || ~isnumeric(Q), error('ssm:ssm_sarimahtd:InputError', 'Q must be a scalar.'); end
if ~isscalar(s) || ~isnumeric(s) || s < 1, error('ssm:ssm_sarimahtd:InputError', 's must be a scalar no less than one.'); end
if nargin < 8, gauss = true; elseif ~islogical(gauss), error('ssm:ssm_sarimahtd:InputError', 'gauss must be logical.'); end
[Z T R Qmat P1 Tmask Rmask Qmask P1mask]    = mat_sarimahtd(p, d, q, P, D, Q, s);
[fun gra psi]                               = fun_sarimahtd(p, d, q, P, D, Q, s, gauss);
if gauss, H = ssmat(0, true); else H = ssdist(1); end
model   = ssmodel(struct('type', 'sarima', 'p', p, 'd', d, 'q', q, 'P', P, 'D', D, 'Q', Q, 's', s), H, Z, ssmat(T, Tmask), ssmat(R, Rmask), ssmat(Qmat, Qmask), ['H' repmat('ng', 1, ~gauss) repmat('T', 1, p+P>0) 'RQP1'], fun, gra, psi, [], ssmat(P1, P1mask));

