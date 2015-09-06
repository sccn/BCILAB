function valid = shouldbevalid(A)

%@SSFUNC/SHOULDBEVALID State space nonlinear function should be valid.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

shouldbevalid(A.ssmat);

f_type          = iscell(A.f) && ndims(A.f) == 2 && (isempty(A.f) || size(A.f, 1) == 1);
df_type         = iscell(A.df) && (isempty(A.df) || isequal(size(A.f), size(A.df)));
horzmask_type   = islogical(A.horzmask) && isequal(size(A.horzmask), [size(A.ssmat, 2) length(A.f)]);
vertmask_type   = islogical(A.vertmask) && isequal(size(A.vertmask), [size(A.ssmat, 1) length(A.f)]);
fmask_type      = islogical(A.fmask) && isequal(size(A.fmask), [1 length(A.f)]);

horzmask        = false(size(A.ssmat, 2), 1);
vertmask        = false(size(A.ssmat, 1), 1);
valid_content   = true;
for i = 1 : length(A.f)
    valid_content   = valid_content ...
        && (isa(A.f{i}, 'function_handle') || (A.fmask(i) && isequal(A.f{i}, 0))) ...
        && (isa(A.df{i}, 'function_handle') || (A.fmask(i) && isequal(A.df{i}, 0))) ...
        && A.ssmat.dmmask(A.vertmask(:, i), A.horzmask(:, i)) ...
        && ~any(horzmask & A.horzmask(:, i)) && ~any(vertmask & A.vertmask(:, i));
    horzmask        = horzmask | A.horzmask(:, i);
    vertmask        = vertmask | A.vertmask(:, i);
end
fmmask                      = false(size(A.ssmat));
fmmask(vertmask, horzmask)  = true;
valid_content               = valid_content && ~any(any(fmmask & A.ssmat.mmask));

m               = size(A.ssmat, 2);
alpha           = randn(m, 1);
t               = floor(100*rand) + 1;
valid_function  = true;
for i = 1 : length(A.f)
    if isa(A.f{i}, 'function_handle')
        y               = A.f{i}(alpha(A.horzmask(:, i)), t);
        valid_function  = valid_function && isequal(size(y), [nnz(A.vertmask(:, i)) 1]);
    end
    if isa(A.df{i}, 'function_handle')
        dy              = A.df{i}(alpha(A.horzmask(:, i)), t);
        valid_function  = valid_function && isequal(size(dy), [nnz(A.horzmask(:, i)) nnz(A.vertmask(:, i))]);
    end
end

type    = f_type && df_type && horzmask_type && vertmask_type && fmask_type;
valid   = type && valid_content && valid_function;

if ~valid, warning('ssm:ssfunc:shouldbevalid', 'Invalid State space nonlinear function.'); end

