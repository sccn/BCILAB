function valid = shouldbevalid(A)

%@SSDIST/SHOULDBEVALID State space non-Gaussian distribution should be valid.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

parent_valid    = shouldbevalid(A.ssmat) && (size(A.ssmat, 1) == size(A.ssmat, 2));

type_type       = isnumeric(A.type) && ndims(A.type) == 2 && (isempty(A.type) || size(A.type, 1) == 1);
matf_type       = iscell(A.matf) && (isempty(A.matf) || isequal(size(A.matf), size(A.type)));
logpf_type      = iscell(A.logpf) && (isempty(A.logpf) || isequal(size(A.logpf), size(A.type)));
diagmask_type   = islogical(A.diagmask) && isequal(size(A.diagmask), [size(A.ssmat, 1) length(A.type)]);
dmask_type      = islogical(A.dmask) && isequal(size(A.dmask), [1 length(A.type)]);

diagmask        = false(size(A.ssmat, 1), 1);
type_content    = true;
matf_content    = true;
logpf_content   = true;
parent_content  = true;
for i = 1 : length(A.type)
    type_content    = type_content && any(A.type(i) == [0 1]);
    matf_content    = matf_content && (isa(A.matf{i}, 'function_handle') || (A.dmask(i) && isequal(A.matf{i}, 0)));
    logpf_content   = logpf_content && (isa(A.logpf{i}, 'function_handle') || (A.dmask(i) && isequal(A.logpf{i}, 0)));
    parent_content  = parent_content && all(all(A.ssmat.dmmask(A.diagmask(:, i), A.diagmask(:, i)))) && ~any(diagmask & A.diagmask(:, i));
    diagmask        = diagmask | A.diagmask(:, i);
end
ngmmask                     = false(size(A.ssmat));
ngmmask(diagmask, diagmask) = true;
parent_content              = parent_content && ~any(any(ngmmask & A.ssmat.mmask));

if A.ssmat.n == 1, n = floor(100*rand + 2); else n = A.ssmat.n; end
p               = size(A.ssmat, 1);
y               = randn(p, n);
theta           = randn(p, n);
eps             = y - theta;
valid_function  = true;
for i = 1 : length(A.type)
    if A.type(i) == 0
        if isa(A.matf{i}, 'function_handle')
            [H ytilde]      = A.matf{i}(abs(round(y(A.diagmask(:, i), :))) + 10, abs(theta(A.diagmask(:, i), :)));
            valid_function  = valid_function && isequal(size(H), [nnz(A.diagmask(:, i))^2 n]) ...
                && isequal(size(ytilde), [nnz(A.diagmask(:, i)) n]);
        end
        if isa(A.logpf{i}, 'function_handle')
            logp            = A.logpf{i}(abs(round(y(A.diagmask(:, i), :))) + 10, abs(theta(A.diagmask(:, i), :)));
            valid_function  = valid_function && isequal(size(logp), [1 n]);
        end
    else % A.type(i) == 1
        if isa(A.matf{i}, 'function_handle')
            H               = A.matf{i}(eps(A.diagmask(:, i), :));
            valid_function  = valid_function && isequal(size(H), [nnz(A.diagmask(:, i))^2 n]);
        end
        if isa(A.logpf{i}, 'function_handle')
            logp            = A.logpf{i}(eps(A.diagmask(:, i), :));
            valid_function  = valid_function && isequal(size(logp), [1 n]);
        end
    end
end

type    = type_type && matf_type && logpf_type && diagmask_type && dmask_type;
content = type_content && matf_content && logpf_content && parent_content;
valid   = parent_valid && type && content && valid_function;

if ~valid, warning('ssm:ssdist:shouldbevalid', 'Invalid state space non-Gaussian distribution.'); end

