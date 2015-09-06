function model = ssm_spline(delta)

%SSM_SPLINE Create SSMODEL object for continuous smoothing spline.
%   model = SSM_SPLINE(delta)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:ssm_spline:NotEnoughInputs', 'Insufficient input arguments.'); end
[Z T Tdmmask Tdvec R]   = mat_spline(delta);
[fun gra psi]           = fun_spline(delta);
model                   = [ssm_gaussian ssmodel('spline', 0, Z, ssmat(T, [], Tdmmask, Tdvec), R, ssmat(zeros(2), [], true(2), zeros(4, 1), true(4, 1)), 'Qd', fun, gra, psi)];
