function res = is_impure_expression(x)
% an impure expression is a MATLAB structure with a .tracking.expression field
res = isfield(x,'tracking') && isfield(x.tracking,'expression');