function res = is_symbol(x)
% check whether an expression is a symbol
res = isequal(exp_head(x),@Symbol);