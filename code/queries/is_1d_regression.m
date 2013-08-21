function x = is_1d_regression(x)
x = isnumeric(x) && size(x,2) == 1;
