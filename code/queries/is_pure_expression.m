function res = is_pure_expression(x)
% a pure expression is anything except for impure expressions
res = ~is_impure_expression(x);