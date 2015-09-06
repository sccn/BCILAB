function model = ssm_reg(x, varname)

%SSM_REG Create SSMODEL object for regression model.
%   model = SSM_REG(x[, varname])
%       x is a m*n matrix storing values of m variables at n time points.
%       varname is optional name for the regression variable.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_reg:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(x), error('ssm:ssm_reg:InputError', 'x must be numeric.'); end
if nargin < 2, varname = 'user'; end
[Z Zdmmask Zdvec T R]   = mat_reg(x, false);
model                   = [ssm_gaussian ssmodel(struct('type', 'regression', 'variable', varname), 0, ssmat(Z, [], Zdmmask, Zdvec), T, R, [])];

