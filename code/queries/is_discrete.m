function x = is_discrete(x)
x = iscell(x) && length(x) == 3 && strcmp(x{1},'disc');
