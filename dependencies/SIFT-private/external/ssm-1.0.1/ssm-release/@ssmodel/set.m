function model = set(model, A, M, func, grad, psi, pmask)

%@SSMODEL/SET Set model elements.
%   model = SET(model0, A, M[, func, grad, psi, pmask])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, error('ssm:ssmodel:set:NotEnoughInputs', 'Insufficient input arguments.'); end
if ischar(A), A = {A}; elseif ~iscell(A), error('ssm:ssmodel:set:InputError', 'A must be a string or a cell array of strings.'); end
if ~isa(M, 'ssmat'), if isnumeric(M), M = ssmat(M); else error('ssm:ssmodel:set:InputError', 'M must be a SSMAT or a matrix.'); end, end
if ~iscell(func), func = {func}; end
if ~iscell(grad), grad = {grad}; end
if nargin < 7, pmask = {true(1, psi.w)}; end

%% Get adjacency masks %%
f           = size(model.A, 1);
adj_M       = logical(eye(18));
adj_allH    = repmat(adj_M(1, :) | adj_M(9, :) | adj_M(15, :), f, 1);
adj_allZ    = repmat(adj_M(2, :) | adj_M(10, :) | adj_M(17, :), f, 1);
adj_allT    = repmat(adj_M(3, :) | adj_M(11, :) | adj_M(18, :), f, 1);
adj_allR    = repmat(adj_M(4, :) | adj_M(12, :), f, 1);
adj_allQ    = repmat(adj_M(5, :) | adj_M(13, :) | adj_M(16, :), f, 1);
adj_allc    = repmat(adj_M(6, :) | adj_M(14, :), f, 1);
adj_a1      = repmat(adj_M(7, :), f, 1);
adj_P1      = repmat(adj_M(8, :), f, 1);

%%%%%%%% TODO: currently only support addition of one element
switch A{1}(1)
    case 'H'
        if ~isequal(size(M), [model.p model.p]) || all(getn(M) ~= [1 model.n]), error('ssm:ssmodel:set:InconsistentModel', 'H must have consistent dimensions.'); end
        model.H                             = M;
        model.A(adj_allH & model.A == 1)    = -1;
    case 'Z'
        if ~isequal(size(M), size(model.Z)) || all(getn(M) ~= [1 model.n]), error('ssm:ssmodel:set:InconsistentModel', 'Z must have consistent dimensions.'); end
        model.Z                             = M;
        model.A(adj_allZ & model.A == 1)    = -1;
    case 'T'
        if ~isequal(size(M), size(model.T)) || all(getn(M) ~= [1 model.n]), error('ssm:ssmodel:set:InconsistentModel', 'T must have consistent dimensions.'); end
        model.T                             = M;
        model.A(adj_allT & model.A == 1)    = -1;
    case 'R'
        if ~isequal(size(M), size(model.R)) || all(getn(M) ~= [1 model.n]), error('ssm:ssmodel:set:InconsistentModel', 'R must have consistent dimensions.'); end
        model.R                             = M;
        model.A(adj_allR & model.A == 1)    = -1;
    case 'Q'
        if ~isequal(size(M), size(model.Q)) || all(getn(M) ~= [1 model.n]), error('ssm:ssmodel:set:InconsistentModel', 'Q must have consistent dimensions.'); end
        model.Q                             = M;
        model.A(adj_allQ & model.A == 1)    = -1;
    case 'c'
        if ~isequal(size(M), size(model.c)) || all(getn(M) ~= [1 model.n]), error('ssm:ssmodel:set:InconsistentModel', 'c must have consistent dimensions.'); end
        model.c                             = M;
        model.A(adj_allc & model.A == 1)    = -1;
    case 'a1'
        if ~isequal(size(M), size(model.a1)) || ~issta(M), error('ssm:ssmodel:set:InconsistentModel', 'a1 must be stationary and have consistent dimensions.'); end
        model.a1                        = M;
        model.A(adj_a1 & model.A == 1)  = -1;
    case 'P1'
        if ~isequal(size(M), size(model.P1)) || ~issta(M), error('ssm:ssmodel:set:InconsistentModel', 'P1 must be stationary and have consistent dimensions.'); end
        model.P1                        = M;
        model.A(adj_P1 & model.A == 1)  = -1;
    otherwise
        error('ssm:ssmodel:set:InputError', 'Unsupported adjacency matrix.');
end

% Remove unused functions and parameters
removal     = all(model.A <= 0, 2);
if any(removal)
    model.A(removal, :)     = [];
    model.func(removal)     = [];
    model.grad(removal)     = [];
    premoval                = sum(vertcat(model.pmask{removal}), 1) > 0; % parameters used by removed functions
    if any(~removal), premoval = premoval & ~(sum(vertcat(model.pmask{~removal}), 1) > 0); end % parameters used by removed functions and not used by any other functions
    model.psi               = remove(model.psi, premoval); % remove said parameters
    model.pmask(removal)    = []; % remove the parameter masks used by removed functions
    f                       = length(model.pmask);
    for j = 1 : f, m.pmask{j}(premoval) = []; end
end

% Add new functions and parameters
model.A             = [model.A; str2adj(A)];
model.func          = [model.func func];
model.grad          = [model.grad grad];
[model.psi pmask2]  = horzcat(model.psi, psi);
model.pmask         = [model.pmask pmask];
for i = 1 : length(model.pmask)
    if i > f, temp = pmask2{2};
    else temp = pmask2{1}; end
    temp(temp)      = model.pmask{i};
    model.pmask{i}  = temp;
end

