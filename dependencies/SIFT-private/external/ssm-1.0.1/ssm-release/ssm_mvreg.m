function model = ssm_mvreg(p, x, dep)

%SSM_MVREG Create SSMODEL object for multivariate regression model.
%   model = SSM_MVREG(p, x[, dep])
%       dep is the p*size(x, 1) dependence logical matrix, each row records
%           the dependence of variable i on regressors in x.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, error('ssm:ssm_mvreg:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(p) || ~isscalar(p), error('ssm:ssm_mvreg:InputError', 'p must be a scalar.'); end
if ~isnumeric(x), error('ssm:ssm_mvreg:InputError', 'x must be a matrix.'); end
if nargin < 3 || isempty(dep), dep = false(0); elseif ~islogical(dep), error('ssm:ssm_mvreg:InputError', 'dep must be a logical matrix.'); end
[Z Zdmmask Zdvec T R]   = mat_mvreg(p, x, dep);
model                   = [ssm_gaussian(p, false) ssmodel(struct('type', 'multivariate regression', 'p', p), zeros(p), ssmat(Z, [], Zdmmask, Zdvec), T, R, [])];

