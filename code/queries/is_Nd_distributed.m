function x = is_Nd_distributed(x)
x = iscell(x) && length(x) == 2 && ischar(x{1}) && strcmp(x{1},'mvn') && iscell(x{2});
