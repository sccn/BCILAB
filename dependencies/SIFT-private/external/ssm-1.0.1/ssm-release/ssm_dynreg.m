function model = ssm_dynreg(x, varname)

%SSM_DYNREG Create SSMODEL object for dynamic regression model.
%   model = SSM_DYNREG(x[, varname])
%       x is a m*n matrix storing values of m variables at n time points.
%       varname is optional name for the regression variable.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_dynreg:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(x), error('ssm:ssm_dynreg:InputError', 'x must be numeric.'); end
if nargin < 3, varname = 'user'; end
[Z Zdmmask Zdvec T R]   = mat_reg(x, true);
[Q Qmmask]              = mat_var(size(Z, 2), false);
[fun gra psi]           = fun_var(size(Z, 2), false, 'xi');
model                   = [ssm_gaussian ssmodel(struct('type', 'dynamic regression', 'variable', varname), 0, ssmat(Z, [], Zdmmask, Zdvec), T, R, ssmat(Q, Qmmask), 'Q', fun, gra, psi)];

