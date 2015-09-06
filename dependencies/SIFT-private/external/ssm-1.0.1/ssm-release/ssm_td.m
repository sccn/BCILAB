function model = ssm_td(y, m, N, td6)

%SSM_TD Create SSMODEL object for trading day component.
%   model = SSM_TD(y, m, N[, td6])
%       y is the starting year.
%       m is the starting month.
%       N is the total number of months.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, error('ssm:ssm_td:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(y) || ~isnumeric(y), error('ssm:ssm_td:InputError', 'y must be a scalar.'); end
if ~isscalar(m) || ~isnumeric(m), error('ssm:ssm_td:InputError', 'm must be a scalar.'); end
if ~isscalar(N) || ~isnumeric(N), error('ssm:ssm_td:InputError', 'N must be a scalar.'); end
if nargin < 4, td6 = true; end
[Y M]                   = ymarray(y, m, N);
[Z Zdmmask Zdvec T R]   = mat_reg(x_td(Y, M, td6), false);
model                   = [ssm_null ssmodel(struct('type', 'trading day', 'nvar', 6, 'y', y, 'm', m, 'N', N), 0, ssmat(Z, [], Zdmmask, Zdvec), T, R, [])];

