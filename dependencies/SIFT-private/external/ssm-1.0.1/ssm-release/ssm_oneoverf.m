function model = ssm_oneoverf(m)

%SSM_ONEOVERF Create SSMODEL object for 1/f noise model.
%   model = SSM_ONEOVERF(m)
%       m is the number of AR terms to use in the approximation.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1 || isempty(m), m = 10; elseif ~isnumeric(m) || ~isscalar(m), error('ssm:ssm_oneoverf:InputError', 'm must be a scalar.'); end
[Z T R P1 Tmmask Rmmask P1mmask]    = mat_arma(m, 0, false);
[fun gra psi]                       = fun_oneoverf(m);
model                               = [ssm_null ssmodel(struct('type', '1/f noise', 'm', m), 0, Z, ssmat(T, Tmmask), ssmat(R, Rmmask), ssmat(0, true), 'TQP1', fun, gra, psi, [], ssmat(P1, P1mmask))];

