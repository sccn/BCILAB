function model = ssm_ee(y, m, N, d)

%SSM_EE Create SSMODEL object for Easter effect component.
%   model = SSM_EE(y, m, N, d)
%       y is the starting year.
%       m is the starting month.
%       N is the total number of months.
%       d is the number of days before Easter.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4, error('ssm:ssm_ee:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isscalar(y) || ~isnumeric(y), error('ssm:ssm_ee:InputError', 'y must be a scalar.'); end
if ~isscalar(m) || ~isnumeric(m), error('ssm:ssm_ee:InputError', 'm must be a scalar.'); end
if ~isscalar(N) || ~isnumeric(N), error('ssm:ssm_ee:InputError', 'N must be a scalar.'); end
if ~isscalar(d) || ~isnumeric(d), error('ssm:ssm_ee:InputError', 'd must be a scalar.'); end
[Y M]                   = ymarray(y, m, N);
[Z Zdmmask Zdvec T R]   = mat_reg(x_ee(Y, M, d), false);
model                   = [ssm_null ssmodel(struct('type', 'easter effect', 'd', d, 'y', y, 'm', m, 'N', N), 0, ssmat(Z, [], Zdmmask, Zdvec), T, R, [])];

