function model = ssm_lom(y, m, N)

%SSM_LOM Create SSMODEL object for length-of-month component.
%   model = SSM_LOM(y, m, N)
%       y is the starting year.
%       m is the starting month.
%       N is the total number of months.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, error('ssm:ssm_lom:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(y) || ~isnumeric(y), error('ssm:ssm_lom:InputError', 'y must be a scalar.'); end
if ~isscalar(m) || ~isnumeric(m), error('ssm:ssm_lom:InputError', 'm must be a scalar.'); end
if ~isscalar(N) || ~isnumeric(N), error('ssm:ssm_lom:InputError', 'N must be a scalar.'); end
[Y M]                   = ymarray(y, m, N);
[Z Zdmmask Zdvec T R]   = mat_reg(eomday(Y, M) - 30.4375, false);
model                   = [ssm_null ssmodel(struct('type', 'length-of-month', 'y', y, 'm', m, 'N', N), 0, ssmat(Z, [], Zdmmask, Zdvec), T, R, [])];

