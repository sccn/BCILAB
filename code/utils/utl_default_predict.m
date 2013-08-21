function y = utl_default_predict(X,M)
% default prediction function, internal to utl_crossval, wraps ml_predict.
% in addition to the formats supported by ml_predict, this function adds support for arrays of the 
% form {NxF,NxD}, the second argument is then interpreted as label vector and omitted
%
% See also:
%   utl_crossval

if iscell(X) && length(X) == 2 && isnumeric(X{1}) && isnumeric(X{2}) && size(X{1},1)==size(X{2},1)
    X = X{1}; end

y = ml_predict(X,M);
