function model = ssmodel(info, H, Z, T, R, Q, AC, func, grad, psi, pmask, P1, a1, c)

%@SSMODEL/SSMODEL State space model class constructor.
%   model = SSMODEL(info, H, Z, T, R, Q) and
%   model = SSMODEL(info, H, Z, T, R, Q, A, func, grad, psi[, pmask,
%       P1, a1, c]) create a general state space model.
%       info is a string, a structure with field 'type' or a cell array of
%           such structures.
%       H is the observation disturbance variance p*p SSMAT or SSDIST.
%       Z is the state to observation p*m SSMAT or SSFUNC.
%       T is the state transition m*m SSMAT or SSFUNC.
%       R is the state disturbance to state m*r SSMAT.
%       Q is the state disturbance variance r*r SSMAT or SSDIST.
%       A is a cell array of strings that specify which part of the model each
%           function updates, each string can be a concatenation (in any
%           order) of the following strings: 'H', 'Z', 'T', 'R', 'Q', 'c',
%           'a1', 'P1', 'Hd', 'Zd', 'Td', 'Rd', 'Qd', 'cd', 'Hng', 'Qng',
%           'Znl', 'Tnl', where suffix 'd', 'ng', 'nl' stands for dynamic,
%           non-Gaussian and nonlinear respectively. Complicated example:
%           ['H' repmat('ng', 1, ~gauss) repmat('T', 1, p+P>0) 'RQP1'].
%       func is a cell array of functions that updates model matrices. If
%           multiple parts of the model are updated the output for each part
%           must be ordered as listed above, but parts not updated can be
%           skipped, for example, a function that updates H, T and Q must
%           output [Hvec Tvec Qvec] in this order.
%       grad is the derivative of func, or [] if not differentiable. For the
%           example function above the gradient function would output [Hvec
%           Tvec Qvec Hgrad Tgrad Qgrad] in this order.
%       psi is a SSPARAM, storing parameter information for the model.
%       pmask is a length(func)*psi.w logical matrix used as parameter masks
%           for the corresponding functions, can be omitted or set to [] if
%           there's only one update function.
%       P1 is the initial state vector variance m*m SSMAT.
%       a1 is the initial state vector mean m*1 SSMAT.
%       c is the state transition constant m*1 SSMAT.
%       All arguments expecting SSMAT can also take numeric matrices.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin == 0
    %% Default constructor %%
    model.name  = '';
    model.Hinfo = [];
    model.info  = {};
    model.mcom  = 0;
    model.H     = ssmat();
    model.Z     = ssmat();
    model.T     = ssmat();
    model.R     = ssmat();
    model.Q     = ssmat();
    model.c     = ssmat(zeros(0, 1));
    model.a1    = ssmat(zeros(0, 1));
    model.P1    = ssmat();
    model.A     = zeros(0, 18);
    model.func  = cell(1, 0);
    model.grad  = cell(1, 0);
    model.psi   = ssparam();
    model.pmask = cell(1, 0);
    model.p     = 0;
    model.n     = 1;
elseif isa(info, 'ssmodel')
    %% Copy constructor %%
    model       = info;
    return;
else
    %%%% Check first group of arguments %%%%
    if nargin < 6, error('ssm:ssmodel:ssmodel:NotEnoughInputs', 'Insufficient input arguments.'); end
    if ischar(info), info = struct('type', info); elseif (~isstruct(info) || ~isfield(info, 'type')) && ~iscell(info), error('ssm:ssmodel:ssmodel:InputError', 'info must be a string, a structure with the field ''type'' or a cell array.'); end
    if ~isa(H, 'ssmat'), if isnumeric(H), H = ssmat(H); else error('ssm:ssmodel:ssmodel:InputError', 'H must be a SSMAT, SSDIST or matrix.'); end, end
    if ~isa(Z, 'ssmat'), if isnumeric(Z), Z = ssmat(Z); else error('ssm:ssmodel:ssmodel:InputError', 'Z must be a SSMAT, SSFUNC or matrix.'); end, end
    if ~isa(T, 'ssmat'), if isnumeric(T), T = ssmat(T); else error('ssm:ssmodel:ssmodel:InputError', 'T must be a SSMAT, SSFUNC or matrix.'); end, end
    if ~isa(R, 'ssmat'), if isnumeric(R), R = ssmat(R); else error('ssm:ssmodel:ssmodel:InputError', 'R must be a SSMAT or matrix.'); end, end
    if ~isa(Q, 'ssmat'), if isnumeric(Q), Q = ssmat(Q); else error('ssm:ssmodel:ssmodel:InputError', 'Q must be a SSMAT, SSDIST or matrix.'); end, end
    p       = size(H, 1);  if size(H, 1) ~= size(H, 2), error('ssm:ssmodel:ssmodel:InconsistentModel', 'H must be square.'); end
    m       = size(Z, 2);  if size(Z, 1) ~= p, error('ssm:ssmodel:ssmodel:InconsistentModel', 'Z must have the same number of rows as H.'); end
    if ~isequal(size(T), [m m]), error('ssm:ssmodel:ssmodel:InconsistentModel', 'T must be square and have the same number of columns as Z.'); end
    r       = size(R, 2);  if size(R, 1) ~= m, error('ssm:ssmodel:ssmodel:InconsistentModel', 'R must have the same number of rows as T.'); end
    if ~isequal(size(Q), [r r]), error('ssm:ssmodel:ssmodel:InconsistentModel', 'Q must be square and have the same number of columns as R.'); end
    %%%%%%%% TODO: Allow for different time durations and pick minimum among > 1's?
    n       = max([H.n Z.n T.n R.n Q.n]); if all(H.n ~= [1 n]) || all(Z.n ~= [1 n]) || all(T.n ~= [1 n]) || all(R.n ~= [1 n]) || all(Q.n ~= [1 n]), error('ssm:ssmodel:ssmodel:InconsistentModel', 'H, Z, T, R, Q must have compatible time durations.'); end

    %%%% Check second group of arguments %%%%
    if nargin >= 7 && ~isempty(AC)
        if nargin < 10, error('ssm:ssmodel:ssmodel:NotEnoughInputs', 'Insufficient input arguments.'); end
        if ischar(AC), AC = {AC}; elseif ~iscell(AC) && ~ischar(AC), error('ssm:ssmodel:ssmodel:InputError', 'A must be a string or cell array of strings.'); end
        if ~iscell(func), error('ssm:ssmodel:ssmodel:InputError', 'func must be a cell array of functions.'); end
        if ~iscell(grad), error('ssm:ssmodel:ssmodel:InputError', 'grad must be a cell array of functions.'); end
        if ~isa(psi, 'ssparam'), error('ssm:ssmodel:ssmodel:InputError', 'psi must be a SSPARAM.'); end
        if nargin >= 11 && ~isempty(pmask), if ~iscell(pmask), error('ssm:ssmodel:ssmodel:InputError', 'pmask must be a cell array of logical row vectors.'); end
        else pmask = {true(1, psi.w)}; end
        f       = length(func); if length(pmask) ~= f, error('ssm:ssmodel:ssmodel:InconsistentModel', 'pmask must have the same length as func.'); end
        if ~isequal(length(AC), f), error('ssm:ssmodel:ssmodel:InconsistentModel', 'A must have the same number of strings as functions.'); end
        if length(grad) ~= f, error('ssm:ssmodel:ssmodel:InconsistentModel', 'grad must have the same length as func.'); end
    else
        AC      = {};
        func    = cell(1, 0);
        grad    = cell(1, 0);
        psi     = ssparam();
        pmask   = cell(1, 0);
    end

    %%%% Check third group of arguments %%%%
    if nargin >= 12 && ~isempty(P1)
        if ~isa(P1, 'ssmat'), if isnumeric(P1), P1 = ssmat(P1); else error('ssm:ssmodel:ssmodel:InputError', 'P1 must be a SSMAT or matrix.'); end, end
        if ~isequal(size(P1), [m m]), error('ssm:ssmodel:ssmodel:InconsistentModel', 'P1 must be a square m*m matrix.'); end
    else P1 = ssmat(diag(repmat(Inf, m, 1)));
    end
    if nargin >= 13 && ~isempty(a1)
        if ~isa(a1, 'ssmat'), if isempty(a1), a1 = ssmat(zeros(m, 1)); elseif isnumeric(a1), a1 = ssmat(a1); else error('ssm:ssmodel:ssmodel:InputError', 'a1 must be a SSMAT or matrix.'); end, end
        if ~isequal(size(a1), [m 1]), error('ssm:ssmodel:ssmodel:InconsistentModel', 'a1 must be a length m column vector.'); end
    else a1 = ssmat(zeros(m, 1));
    end
    if nargin >= 14
        if ~isa(c, 'ssmat'), if isempty(c), c = ssmat(zeros(m, 1)); elseif isnumeric(c), c = ssmat(c); else error('ssm:ssmodel:ssmodel:InputError', 'c must be a SSMAT or matrix.'); end, end
        if ~isequal(size(c), [m 1]), error('ssm:ssmodel:ssmodel:InconsistentModel', 'c must be a length m column vector.'); end
        if c.n ~= [1 n], error('ssm:ssmodel:ssmodel:InconsistentModel', 'c must have compatible time duration.'); end
    else c = ssmat(zeros(m, 1));
    end

    %%%% Construct model %%%%
    model.name  = '';
    if m > 0
        model.Hinfo = [];
        if ~iscell(info), model.info = {info};
        else model.info = info; end
        model.mcom  = [0 m];
    else
        model.Hinfo = info;
        model.info  = cell(1, 0);
        model.mcom  = 0;
    end
    model.H     = H;
    model.Z     = Z;
    model.T     = T;
    model.R     = R;
    model.Q     = Q;
    model.c     = c;
    model.a1    = a1;
    model.P1    = P1;
    model.A     = str2adj(AC);
    model.func  = func;
    model.grad  = grad;
    model.psi   = psi;
    model.pmask = pmask;
    model.p     = p;
    model.n     = n;
end

%% Register object instance %%
model   = class(model, 'ssmodel');

