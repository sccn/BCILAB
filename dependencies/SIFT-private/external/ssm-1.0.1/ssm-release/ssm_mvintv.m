function model = ssm_mvintv(p, n, type, tau)

%SSM_MVINTV Create SSMODEL object for multivariate intervention component.
%   model = SSM_MVINTV(p, n, type, tau)
%       type can be 'step', 'pulse', 'slope' or 'null' or a cell array.
%       tau is the onset time or an array of onset times.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4, error('ssm:ssm_mvintv:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(p) || ~isscalar(p), error('ssm:ssm_mvintv:InputError', 'p must be a scalar.'); end
if ~isnumeric(n) || ~isscalar(n), error('ssm:ssm_mvintv:InputError', 'n must be a scalar.'); end
if ~ischar(type) && ~iscell(type), error('ssm:ssm_mvintv:InputError', 'type must be a string or a cell array of strings.'); end
if ~isnumeric(tau) || (~isscalar(tau) && length(tau)~=p), error('ssm:ssm_mvintv:InputError', 'tau must be a scalar or a p vector.'); end
if ischar(type), [temp{1:p}] = deal(type); type = temp; end
if isscalar(tau), tau = repmat(tau, 1, p); end
[Z Zdmmask Zdvec T R]   = mat_mvintv(p, n, type, tau);
model                   = [ssm_null(p) ssmodel(struct('type', 'multivariate intervention', 'p', p, 'subtype', {type}, 'tau', tau), zeros(p), ssmat(Z, [], Zdmmask, Zdvec), T, R, [])];

