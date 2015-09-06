function model = ssm_lpt(d)

%SSM_LPT Create SSMODEL object for local polynomial trend model.
%   model = SSM_LPT(d)
%       d is the order of the polynomial trend.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_lpt:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(d) || ~isnumeric(d), error('ssm:ssm_lpt:InputError', 'd must be a scalar.'); end
[Z T R]         = mat_lpt(d, true);
[Q Qmmask]      = mat_var(d + 1, false);
[fun gra psi]   = fun_var(d + 1, false, 'zeta');
model           = [ssm_gaussian ssmodel(struct('type', 'local polynomial trend', 'd', d), 0, Z, T, R, ssmat(Q, Qmmask), 'Q', fun, gra, psi)];

