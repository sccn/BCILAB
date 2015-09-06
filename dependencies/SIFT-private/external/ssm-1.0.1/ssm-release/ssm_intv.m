function model = ssm_intv(n, type, tau)

%SSM_INTV Create SSMODEL object for intervention component.
%   model = SSM_INTV(n, type, tau)
%       n is the total time duration.
%       type can be 'step', 'pulse', 'slope' or 'null'.
%       tau is the onset time.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, error('ssm:ssm_intv:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(n) || ~isscalar(n), error('ssm:ssm_intv:InputError', 'n must be a scalar.'); end
if ~ischar(type), error('ssm:ssm_intv:InputError', 'type must be a string.'); end
if ~isnumeric(tau) || ~isscalar(tau), error('ssm:ssm_intv:InputError', 'tau must be a scalar.'); end
[Z Zdmmask Zdvec T R]   = mat_reg(x_intv(n, type, tau), false);
model                   = [ssm_null ssmodel(struct('type', 'intervention', 'subtype', type, 'tau', tau), 0, ssmat(Z, [], Zdmmask, Zdvec), T, R, [])];

