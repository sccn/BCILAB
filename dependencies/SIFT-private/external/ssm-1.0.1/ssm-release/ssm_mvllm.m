function model = ssm_mvllm(p, cov)

%SSM_MVLLM Create SSMODEL object for multivariate local level model.
%   model = SSM_MVLLM(p[, cov])
%       p is the number of variables.
%       cov specifies complete covariance if true, or complete independence if
%           false, extended to a vector where needed.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_mvllm:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(p) || ~isscalar(p), error('ssm:ssm_mvllm:InputError', 'p must be a scalar.'); end
if nargin < 2 || isempty(cov), cov = true(1, 2); elseif ~islogical(cov), error('ssm:ssm_mvllm:InputError', 'cov must be logical.'); elseif isscalar(cov), cov = [cov cov]; end
[Q Qmmask]      = mat_var(p, cov(2));
[fun gra psi]   = fun_var(p, cov(2), 'zeta');
model           = [ssm_gaussian(p, cov(1)) ssmodel(struct('type', 'multivariate local level', 'p', p), zeros(p), eye(p), eye(p), eye(p), ssmat(Q, Qmmask), 'Q', fun, gra, psi)];

