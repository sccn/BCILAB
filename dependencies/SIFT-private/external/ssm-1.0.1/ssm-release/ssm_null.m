function model = ssm_null(p)

%SSM_NULL Create SSMODEL object for null/no noise model.
%   model = SSM_NULL(p)
%       p is the number of variables.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, p = 1; elseif ~isnumeric(p) || ~isscalar(p), error('ssm:ssm_null:InputError', 'p must be a scalar.'); end
model   = ssmodel(struct('type', 'null noise', 'p', p), zeros(p), zeros(p, 0), [], [], []);

