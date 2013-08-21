function res = has_function_symbol(x)
% internal: a function_handle from which we can obtain a function symbol: every function handle except for non-symbolic lambdas.
res = isa(x,'function_handle') && (isempty(strmatch('@',char(x))) || is_symbolic_lambda(x));
