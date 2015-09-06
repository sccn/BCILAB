function model = ssm_irw(d)

%SSM_IRW Create SSMODEL object for integrated random walk model.
%   model = SSM_IRW(d)
%       d is the order of integration.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_irw:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(d) || ~isnumeric(d), error('ssm:ssm_irw:InputError', 'd must be a scalar.'); end
[Z T R]         = mat_lpt(d, false);
[Q Qmmask]      = mat_var(1, false);
[fun gra psi]   = fun_var(1, false, 'zeta');
model           = [ssm_gaussian ssmodel(struct('type', 'integrated random walk', 'd', d), 0, Z, T, R, ssmat(Q, Qmmask), 'Q', fun, gra, psi)];

