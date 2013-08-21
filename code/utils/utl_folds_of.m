function result = utl_folds_of(scheme)
% Get the number of folds for a given cross-validation scheme.
% For 'loo', this is assumed to be ~10.
%
% See also:
%   utl_crossval

if isnumeric(scheme)    
    result = prod(scheme(1:min(2,end)));
elseif isequal(scheme,{})
    % this is the guessing part
    result = 10;
elseif strcmp(scheme,'loo')
    % same here
    result = 10;
elseif iscell(scheme) && any(strcmp(scheme{1},{'chron','block'}))
    result = scheme{2};
elseif iscell(scheme) && any(strcmp(scheme{1},{'subchron','subblock'}))
    result = scheme{3};
elseif strcmp(scheme,'trainerr')
    result = 1;
else
    error('Unrecognized cross-validation scheme format.');
end
