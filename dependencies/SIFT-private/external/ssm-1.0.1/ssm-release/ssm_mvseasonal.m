function model = ssm_mvseasonal(p, cov, type, s)

%SSM_MVSEASONAL Create SSMODEL object for multivariate seasonal component.
%   model = SSM_MVSEASONAL(p, cov, type, s)
%       p is the number of variables.
%       cov specifies complete covariance if true, or complete independence if
%           false, extended to a vector where needed.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4, error('ssm:ssm_mvseasonal:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(p) || ~isscalar(p), error('ssm:ssm_mvseasonal:InputError', 'p must be a scalar.'); end
if isempty(cov), cov = true; elseif ~islogical(cov), error('ssm:ssm_mvseasonal:InputError', 'cov must be logical.'); end
if ~ischar(type), error('ssm:ssm_mvseasonal:InputError', 'type must be a string.'); end
if ~isnumeric(s) || ~isscalar(s), error('ssm:ssm_mvseasonal:InputError', 's must be a scalar.'); end
switch type
    case 'dummy'
        [Z T R]         = mat_mvdummy(p, s, false);
        [Q Qmmask]      = mat_var(p, cov);
        [fun gra psi]   = fun_var(p, cov, 'omega');
    case 'dummy fixed'
        [Z T R]         = mat_mvdummy(p, s, true);
        Q               = [];
        fun             = {};
    case 'h&s'
        % Multivariate H&S seasonal is always assumed independent
        [Z T R]         = mat_mvhs(p, s);
        [Q Qmmask]      = mat_wvar(p, s);
        [fun gra psi]   = fun_wvar(p, s, 'omega');
    case 'trig1'
        [Z T R]         = mat_mvtrig(p, s, false);
        [Q Qmmask]      = mat_dupvar(p, cov, s - 1);
        [fun gra psi]   = fun_dupvar(p, cov, s - 1, 'omega');
    case 'trig2'
        [Z T R]         = mat_mvtrig(p, s, false);
        [Q Qmmask]      = mat_interlvar(p, s - 1, cov);
        [fun gra psi]   = fun_interlvar(p, s - 1, cov, 'omega');
    case 'trig fixed'
        [Z T R]         = mat_mvtrig(p, s, true);
        Q               = [];
        fun             = {};
    otherwise
        error('ssm:ssm_mvseasonal:NotSupported', 'Unsupported multivariate seasonal type.');
end
if isempty(fun), model = [ssm_null(p) ssmodel(struct('type', 'multivariate seasonal', 'p', p, 'subtype', type, 's', s), zeros(p), Z, T, R, Q)];
else model = [ssm_null(p) ssmodel(struct('type', 'multivariate seasonal', 'p', p, 'subtype', type, 's', s), zeros(p), Z, T, R, ssmat(Q, Qmmask), 'Q', fun, gra, psi)];
end

