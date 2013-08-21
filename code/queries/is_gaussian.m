function x = is_gaussian(x)
x = iscell(x) && length(x) == 3 && strcmp(x{1},'norm');
