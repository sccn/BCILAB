function x = is_custom_distributed(x)
x = iscell(x) && length(x) == 2 && ischar(x{1}) && isa(x{1},'function_handle') && iscell(x{2});
