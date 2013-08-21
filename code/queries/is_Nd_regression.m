function x = is_Nd_regression(x)
x = isnumeric(x) && size(x,2) > 1;
