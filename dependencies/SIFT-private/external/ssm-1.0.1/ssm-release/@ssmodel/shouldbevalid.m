function valid = shouldbevalid(A)

%@SSMODEL/SHOULDBEVALID State space model should be valid.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

shouldbevalid(A.H);
shouldbevalid(A.Z);
shouldbevalid(A.T);
shouldbevalid(A.R);
shouldbevalid(A.Q);
shouldbevalid(A.c);
shouldbevalid(A.a1);
shouldbevalid(A.P1);

valid_type  = ischar(A.name) && (isempty(A.Hinfo) || (isstruct(A.Hinfo) && isfield(A.Hinfo, 'type'))) ...
    && iscell(A.info) && ndims(A.info) == 2 && size(A.info, 1) == 1 ...
    && isnumeric(A.mcom) && ndims(A.mcom) == 2 && size(A.mcom, 1) == 1 ...
    && isa(A.H, 'ssmat') && isa(A.Z, 'ssmat') && isa(A.T, 'ssmat') ...
    && isa(A.R, 'ssmat') && isa(A.Q, 'ssmat') && isa(A.c, 'ssmat') ...
    && isa(A.a1, 'ssmat') && isa(A.P1, 'ssmat') ...
    && isnumeric(A.A) && ndims(A.A) == 2 && size(A.A, 2) == 18 ...
    && iscell(A.func) && iscell(A.grad) && isa(A.psi, 'ssparam') ...
    && iscell(A.pmask) && isnumeric(A.p) && isscalar(A.p) && isnumeric(A.n) && isscalar(A.n);
valid_dim   = size(A.H, 1) == size(A.H, 2) && size(A.H, 2) == size(A.Z, 1) ...
    && size(A.Z, 2) == size(A.T, 1) && size(A.T, 1) == size(A.T, 2) && size(A.T, 2) == size(A.R, 1) ...
    && size(A.R, 2) == size(A.Q, 1) && size(A.Q, 1) == size(A.Q, 2);

p   = size(A.H, 1);
m   = size(A.T, 1);
n   = max([A.H.n A.Z.n A.T.n A.R.n A.Q.n A.c.n]);

valid_dim   = valid_dim && isequal(size(A.c), [m 1]) && isequal(size(A.a1), [m 1]) && isequal(size(A.P1), [m m]) ...
    && any(A.H.n == [1 n]) && any(A.Z.n == [1 n]) && any(A.T.n == [1 n]) && any(A.R.n == [1 n]) ...
    && any(A.Q.n == [1 n]) && any(A.c.n == [1 n]) && A.a1.n == 1 && A.P1.n == 1 && A.p == p && A.n == n;

mcom            = diff(A.mcom);
valid_content   = all(mcom > 0) && sum(mcom) == m && length(A.info) == length(mcom);
for i = 1 : length(A.info)
    valid_content   = valid_content && isstruct(A.info{i}) && isfield(A.info{i}, 'type');
end
f   = size(A.func, 2);
w   = A.psi.w;
valid_content   = valid_content && size(A.A, 1) == f && ndims(A.func) == 2 && size(A.func, 1) == 1 ...
    && isequal(size(A.func), size(A.grad)) && isequal(size(A.func), size(A.pmask));

X               = 5*randn(1, w);
valid_function  = true;
valid_content2  = true;
Result          = cell(f, 18);
for i = 1 : f
    valid_content2  = valid_content2 && isa(A.func{i}, 'function_handle') ...
        && (isa(A.grad{i}, 'function_handle') || isempty(A.grad{i})) ...
        && islogical(A.pmask{i}) && isequal(size(A.pmask{i}), [1 w]);
    [Result{i, A.A(i, :) ~= 0}]   = A.func{i}(X(A.pmask{i}));
    if ~isempty(A.grad{i})
        [ResultF{1:nnz(A.A(i, :) ~= 0)}]    = Result{i, A.A(i, :) ~= 0};
        [ResultF2{1:nnz(A.A(i, :) ~= 0)} ResultG{1:nnz(A.A(i, :) ~= 0)}]    = A.grad{i}(X(A.pmask{i}));
        for j = 1 : length(ResultF)
            valid_function  = valid_function && isequal(ResultF, ResultF2) ...
                && size(ResultF{j}, 1) == size(ResultG{j}, 1) && size(ResultG{j}, 2) == nnz(A.pmask{i});
        end
    end
end
for i = 1 : 18, Result{1, i} = vertcat(Result{:, i}); end

%% Adjacency logical matrix entries %%
adj_M   = logical(eye(18));
adj_H   = adj_M(1, :);
adj_Z   = adj_M(2, :);
adj_T   = adj_M(3, :);
adj_R   = adj_M(4, :);
adj_Q   = adj_M(5, :);
adj_c   = adj_M(6, :);
adj_a1  = adj_M(7, :);
adj_P1  = adj_M(8, :);
adj_Hd  = adj_M(9, :);
adj_Zd  = adj_M(10, :);
adj_Td  = adj_M(11, :);
adj_Rd  = adj_M(12, :);
adj_Qd  = adj_M(13, :);
adj_cd  = adj_M(14, :);
adj_Hng = adj_M(15, :);
adj_Qng = adj_M(16, :);
adj_Znl = adj_M(17, :);
adj_Tnl = adj_M(18, :);

if any(any(A.H.mmask)), valid_function = valid_function && isnumeric(Result{1, adj_H}) && isequal(size(Result{1, adj_H}), [nnz(A.H.mmask) 1]); end
if any(any(A.Z.mmask)), valid_function = valid_function && isnumeric(Result{1, adj_Z}) && isequal(size(Result{1, adj_Z}), [nnz(A.Z.mmask) 1]); end
if any(any(A.T.mmask)), valid_function = valid_function && isnumeric(Result{1, adj_T}) && isequal(size(Result{1, adj_T}), [nnz(A.T.mmask) 1]); end
if any(any(A.R.mmask)), valid_function = valid_function && isnumeric(Result{1, adj_R}) && isequal(size(Result{1, adj_R}), [nnz(A.R.mmask) 1]); end
if any(any(A.Q.mmask)), valid_function = valid_function && isnumeric(Result{1, adj_Q}) && isequal(size(Result{1, adj_Q}), [nnz(A.Q.mmask) 1]); end
if any(any(A.c.mmask)), valid_function = valid_function && isnumeric(Result{1, adj_c}) && isequal(size(Result{1, adj_c}), [nnz(A.c.mmask) 1]); end
if any(any(A.a1.mmask)), valid_function = valid_function && isnumeric(Result{1, adj_a1}) && isequal(size(Result{1, adj_a1}), [nnz(A.a1.mmask) 1]); end
if any(any(A.P1.mmask)), valid_function = valid_function && isnumeric(Result{1, adj_P1}) && isequal(size(Result{1, adj_P1}), [nnz(A.P1.mmask) 1]); end
if any(any(A.H.dvmask)), valid_function = valid_function && isnumeric(Result{1, adj_Hd}) && size(Result{1, adj_Hd}, 1) == nnz(A.H.dvmask); end
if any(any(A.Z.dvmask)), valid_function = valid_function && isnumeric(Result{1, adj_Zd}) && size(Result{1, adj_Zd}, 1) == nnz(A.Z.dvmask); end
if any(any(A.T.dvmask)), valid_function = valid_function && isnumeric(Result{1, adj_Td}) && size(Result{1, adj_Td}, 1) == nnz(A.T.dvmask); end
if any(any(A.R.dvmask)), valid_function = valid_function && isnumeric(Result{1, adj_Rd}) && size(Result{1, adj_Rd}, 1) == nnz(A.R.dvmask); end
if any(any(A.Q.dvmask)), valid_function = valid_function && isnumeric(Result{1, adj_Qd}) && size(Result{1, adj_Qd}, 1) == nnz(A.Q.dvmask); end
if any(any(A.c.dvmask)), valid_function = valid_function && isnumeric(Result{1, adj_cd}) && size(Result{1, adj_cd}, 1) == nnz(A.c.dvmask); end
if isa(A.H, 'ssdist') && any(A.H.dmask), valid_function = valid_function && iscell(Result{1, adj_Hng}) && isequal(size(Result{1, adj_Hng}), [nnz(A.H.dmask) 2]); end
if isa(A.Q, 'ssdist') && any(A.Q.dmask), valid_function = valid_function && iscell(Result{1, adj_Qng}) && isequal(size(Result{1, adj_Qng}), [nnz(A.Q.dmask) 2]); end
if isa(A.Z, 'ssfunc') && any(A.Z.fmask), valid_function = valid_function && iscell(Result{1, adj_Znl}) && isequal(size(Result{1, adj_Znl}), [nnz(A.Z.fmask) 2]); end
if isa(A.T, 'ssfunc') && any(A.T.fmask), valid_function = valid_function && iscell(Result{1, adj_Tnl}) && isequal(size(Result{1, adj_Tnl}), [nnz(A.T.fmask) 2]); end

valid   = valid_type && valid_dim && valid_content && valid_content2 && valid_function;

if ~valid, warning('ssm:ssmodel:shouldbevalid', 'Invalid state space model.'); end

