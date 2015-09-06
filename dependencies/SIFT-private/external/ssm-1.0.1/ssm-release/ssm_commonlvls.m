function model = ssm_commonlvls(p, A_ast, a_ast)

%SSM_COMMONLVLS Create SSMODEL object for common levels model.
%   model = SSM_COMMONLVLS(p, A_ast, a_ast)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, error('ssm:ssm_commonlvls:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(p) || ~isscalar(p), error('ssm:ssm_commonlvls:InputError', 'p must be a scalar.'); end
if ~isnumeric(A_ast), error('ssm:ssm_commonlvls:InputError', 'A_ast must be a matrix.'); end
if ~isnumeric(a_ast), error('ssm:ssm_commonlvls:InputError', 'a_ast must be a matrix.'); end
[Z T R a1 P1]   = mat_commonlvls(p, A_ast, a_ast);
[Q Qmmask]      = mat_var(size(R, 2), true);
[fun gra psi]   = fun_var(size(R, 2), true, 'zeta');
model           = [ssm_gaussian(p, true) ssmodel('common levels', zeros(p), Z, T, R, ssmat(Q, Qmmask), 'Q', fun, gra, psi, [], P1, a1)];

