function res = is_dataset(x)
% Determine whether some object is a data set.
res = all(isfield(x,{'data','srate'}));

