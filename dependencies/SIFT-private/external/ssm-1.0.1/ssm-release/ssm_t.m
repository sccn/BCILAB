function model = ssm_t(nu)

%SSM_T Create SSMODEL object for t-distribution noise model.
%   model = SSM_T([nu])
%       nu is the degree of freedom.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, nu = 0; elseif ~isnumeric(nu) || ~isscalar(nu), error('ssm:ssm_t:InputError', 'nu must be a scalar.'); end
[fun grad param]    = fun_t(nu);
model               = ssmodel(struct('type', 't-distribution noise', 'nu', nu), ssdist(1), zeros(1, 0), [], [], [], 'Hng', fun, grad, param);

