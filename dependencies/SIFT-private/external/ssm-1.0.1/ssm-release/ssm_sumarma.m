function model = ssm_sumarma(p, q, D, s, mean)

%SSM_SUMARMA Create SSMODEL object for sum integrated ARMA(p, q) model of order s.
%   model = SSM_SUMARMA(p, q, D, s[, mean])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4, error('ssm:ssm_sumarma:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(p) || ~isnumeric(p), error('ssm:ssm_sumarma:InputError', 'p must be a scalar.'); end
if ~isscalar(q) || ~isnumeric(q), error('ssm:ssm_sumarma:InputError', 'q must be a scalar.'); end
if ~isscalar(D) || ~isnumeric(D), error('ssm:ssm_sumarma:InputError', 'D must be a scalar.'); end
if ~isscalar(s) || ~isnumeric(s), error('ssm:ssm_sumarma:InputError', 's must be a scalar.'); end
if nargin < 5, mean = false; end
[Z T R P1 Tmask Rmask P1mask]   = mat_sumarma(p, q, D, s, mean);
[fun gra psi]                   = fun_arma(p, q);
model                           = [ssm_null ssmodel(struct('type', 'sumarma', 'p', p, 'q', q, 'D', D, 's', s), 0, Z, ssmat(T, Tmask), ssmat(R, Rmask), ssmat(0, true), [repmat('T', 1, p>0) repmat('R', 1, q>0) 'QP1'], fun, gra, psi, [], ssmat(P1, P1mask))];

