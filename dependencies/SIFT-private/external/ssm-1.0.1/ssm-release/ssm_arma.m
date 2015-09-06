function model = ssm_arma(p, q, mean)

%SSM_ARMA Create SSMODEL object for ARMA(p, q) model.
%   model = SSM_ARMA(p, q[, mean])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, error('ssm:ssm_arma:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(p) || ~isnumeric(p), error('ssm:ssm_arma:InputError', 'p must be a scalar.'); end
if ~isscalar(q) || ~isnumeric(q), error('ssm:ssm_arma:InputError', 'q must be a scalar.'); end
if nargin < 3, mean = false; end
[Z T R P1 Tmask Rmask P1mask]   = mat_arma(p, q, mean);
[fun gra psi]                   = fun_arma(p, q);
model                           = [ssm_null ssmodel(struct('type', 'arima', 'p', p, 'd', 0, 'q', q), 0, Z, ssmat(T, Tmask), ssmat(R, Rmask), ssmat(0, true), [repmat('T', 1, p>0) repmat('R', 1, q>0) 'QP1'], fun, gra, psi, [], ssmat(P1, P1mask))];

