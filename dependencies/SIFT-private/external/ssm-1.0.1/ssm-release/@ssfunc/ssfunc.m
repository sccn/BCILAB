function func = ssfunc(p, m, f, df)

%@SSFUNC/SSFUNC State space nonlinear function class constructor.
%   func = SSFUNC(p, m[, f, df])
%       f is a function that maps m*1 vectors and time t to p*1 vectors.
%       df is the derivative of f that maps m*1 vectors and time t to p*m
%           matrices.
%       (if f and df are not both provided, the nonlinear function is assumed
%           to be variable)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin == 0
    %% Default constructor %%
    parent          = ssmat;
    func.f          = {};
    func.df         = {};
    func.horzmask   = false(0);
    func.vertmask   = false(0);
    func.fmask      = false(1, 0);
elseif isa(p, 'ssfunc')
    %% Copy constructor %%
    func            = p;
    return;
else
    if isa(p, 'ssmat')
        parent          = p;
        func.f          = {};
        func.df         = {};
        func.horzmask   = false(size(p, 2), 0);
        func.vertmask   = false(size(p, 1), 0);
        func.fmask      = false(1, 0);
    else
        if nargin < 2, error('ssm:ssfunc:ssfunc:NotEnoughInputs', 'Insufficient input arguments.'); end
        if ~isnumeric(p) || ~isscalar(p), error('ssm:ssfunc:ssfunc:InputError', 'p must be a scalar.'); end
        if ~isnumeric(m) || ~isscalar(m), error('ssm:ssfunc:ssfunc:InputError', 'm must be a scalar.'); end
        if nargin < 4
            %% Variable general nonlinear function %%
            parent          = ssmat(zeros(p, m), [], true(p, m), zeros(p*m, 1));
            func.f          = {0};
            func.df         = {0};
            func.horzmask   = true(m, 1);
            func.vertmask   = true(p, 1);
            func.fmask      = true;
        else
            %% Constant general nonlinear function %%
            if ~isa(f, 'function_handle') && ~isequal(f, 0), error('ssm:ssfunc:ssfunc:InputError', 'f must be a function.'); end
            if ~isa(df, 'function_handle') && ~isequal(df, 0), error('ssm:ssfunc:ssfunc:InputError', 'df must be a function.'); end
            parent          = ssmat(zeros(p, m), [], true(p, m), zeros(p*m, 1));
            func.f          = {f};
            func.df         = {df};
            func.horzmask   = true(m, 1);
            func.vertmask   = true(p, 1);
            func.fmask      = false;
        end
    end
end

superiorto('ssmat');

%% Register object instance %%
func    = class(func, 'ssfunc', parent);

