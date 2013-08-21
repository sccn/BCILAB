function x = is_1d_distributed(x)
x = iscell(x) && length(x) == 2 && ischar(x{1}) && isnumeric(x{2});
