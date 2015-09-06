function model = ssm_gaussian(p, cov)

%SSM_GAUSSIAN Create SSMODEL object for Gaussian noise model.
%   model = SSM_GAUSSIAN(p[, cov])
%       p is the number of variables.
%       cov specifies complete covariance if true, or complete independence if
%           false.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, p = 1; elseif ~isnumeric(p) || ~isscalar(p) || p < 1, error('ssm:ssm_gaussian:InputError', 'p must be a scalar larger than 1.'); end
if nargin < 2, cov = true; elseif ~islogical(cov) || ~isscalar(cov), error('ssm:ssm_gaussian:InputError', 'cov must be a logical scalar.'); end
[H Hmmask]      = mat_var(p, cov);
[fun gra psi]   = fun_var(p, cov, 'epsilon');
model           = ssmodel(struct('type', 'Gaussian noise', 'p', p), ssmat(H, Hmmask), zeros(p, 0), [], [], [], 'H', fun, gra, psi);

