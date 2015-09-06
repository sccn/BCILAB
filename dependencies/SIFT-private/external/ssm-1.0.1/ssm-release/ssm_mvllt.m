function model = ssm_mvllt(p, cov)

%SSM_MVLLT Create SSMODEL object for multivariate local level trend model.
%   model = SSM_MVLLT(p[, cov])
%       p is the number of variables.
%       cov specifies complete covariance if true, or complete independence if
%           false, extended to a vector where needed.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_mvllt:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(p) || ~isscalar(p), error('ssm:ssm_mvllt:InputError', 'p must be a scalar.'); end
if nargin < 2 || isempty(cov), cov = true(1, 3); elseif ~islogical(cov), error('ssm:ssm_mvllt:InputError', 'cov must be logical.'); elseif isscalar(cov), cov = [cov cov cov]; end
[Q Qmmask]      = mat_interlvar(p, 2, cov(2:3));
[fun gra psi]   = fun_interlvar(p, 2, cov(2:3), {'xi' 'zeta'});
[Z T R]         = mat_mvllt(p);
model           = [ssm_gaussian(p, cov(1)) ssmodel(struct('type', 'multivariate local linear trend', 'p', p), zeros(p), Z, T, R, ssmat(Q, Qmmask), 'Q', fun, gra, psi)];

