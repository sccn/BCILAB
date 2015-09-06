function model = ssm_mvcycle(p, cov)

%SSM_MVCYCLE Create SSMODEL object for multivariate cycle component.
%   model = SSM_MVCYCLE(p[, cov])
%       p is the number of variables.
%       cov specifies complete covariance if true, or complete independence if
%           false, extended to a vector where needed.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_mvcycle:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(p) || ~isscalar(p), error('ssm:ssm_mvcycle:InputError', 'p must be a scalar.'); end
if nargin < 2 || isempty(cov), cov = true; elseif ~islogical(cov), error('ssm:ssm_mvcycle:InputError', 'cov must be logical.'); end
[Z T Tmmask R]      = mat_mvcycle(p);
[Tfun Tgra psi{1}]  = fun_mvcycle(p);
[Q Qmmask]          = mat_dupvar(p, cov, 2);
[Qfun Qgra psi{2}]  = fun_dupvar(p, cov, 2, 'omega tilde');
[psi pmask]         = horzcat(psi{:});
model               = [ssm_null(p) ssmodel('multivariate sine cycle', zeros(p), Z, ssmat(T, Tmmask), R, ssmat(Q, Qmmask), {'T' 'Q'}, [Tfun Qfun], [Tgra Qgra], psi, pmask)];

