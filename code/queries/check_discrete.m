function x = check_discrete(x)
x = is_discrete(x) && ischar(x{1}) && strcmp(x{1},'disc') && isnumeric(x{2}) && isnumeric(x{3}) && size(x{2},2) == size(x{3},1);
