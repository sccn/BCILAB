function model = ssm_seasonal(type, s)

%SSM_SEASONAL Create SSMODEL object for seasonal components.
%   model = SSM_SEASONAL(type, s)
%       type can be 'dummy', 'dummy fixed', 'h&s', 'trig1', 'trig2' or 'trig
%           fixed'.
%       s is the seasonal period.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, error('ssm:ssm_seasonal:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~ischar(type), error('ssm:ssm_seasonal:InputError', 'type must be a string.'); end
if ~isnumeric(s) || ~isscalar(s), error('ssm:ssm_seasonal:InputError', 's must be a scalar.'); end
switch type
    case 'dummy' %% The dummy seasonal component %%
        [Z T R]         = mat_dummy(s, false);
        [Q Qmmask]      = mat_var(1, false);
        [fun gra psi]   = fun_var(1, false, 'omega');
    case 'dummy fixed' %% The constant dummy seasonal component %%
        [Z T R]         = mat_dummy(s, true);
        Q               = [];
        fun             = {};
    case 'h&s' %% The seasonal component suggested by Harrison and Stevens (1976) %%
        [Z T R]         = mat_hs(s);
        [Q Qmmask]      = mat_wvar(1, s);
        [fun gra psi]   = fun_wvar(1, s, 'omega');
    case 'trig1' %% Trigonometric seasonal component with one variance %%
        [Z T R]         = mat_trig(s, false);
        [Q Qmmask]      = mat_dupvar(1, false, s - 1);
        [fun gra psi]   = fun_dupvar(1, false, s - 1, 'omega');
    case 'trig2' %% Trigonometric seasonal component with full variance %%
        [Z T R]         = mat_trig(s, false);
        [Q Qmmask]      = mat_var(s-1, false);
        [fun gra psi]   = fun_var(s-1, false, 'omega');
    case 'trig fixed' %% Constant trigonometric seasonal component %%
        [Z T R]         = mat_trig(s, true);
        Q               = [];
        fun             = {};
    otherwise
        error('ssm:ssm_seasonal:NotSupported', 'Unsupported seasonal type.');
end
if isempty(fun), model = [ssm_null ssmodel(struct('type', 'seasonal', 'subtype', type, 's', s), 0, Z, T, R, Q)];
else model = [ssm_null ssmodel(struct('type', 'seasonal', 'subtype', type, 's', s), 0, Z, T, R, ssmat(Q, Qmmask), 'Q', fun, gra, psi)];
end

