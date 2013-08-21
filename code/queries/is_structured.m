function x = is_structured(x)
x = iscell(x) && length(x) == 2 && ischar(x{1}) && strcmp(x{1},'struct') && iscell(x{2});
