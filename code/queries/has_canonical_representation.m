function res = has_canonical_representation(x)
% determine whether an expression is represented as a struct with the fields 'head' and 'parts'.
res = all(isfield(x,{'head','parts'}));
