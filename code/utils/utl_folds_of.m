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
    if length(scheme)<2 || isempty(scheme{2})
        error('The given k-fold chron/block cross-validation scheme is lacking the value for k (second cell element).'); end
    result = scheme{2};
elseif iscell(scheme) && any(strcmp(scheme{1},{'subchron','subblock'}))
    if length(scheme)<3 || isempty(scheme{3})
        error('The given k-fold subchron/subblock cross-validation scheme is lacking the value for k (third cell element).'); end
    result = scheme{3};
elseif strcmp(scheme,'trainerr')
    result = 1;
else
    error('Unrecognized cross-validation scheme format (see utl_crossval for allowed formats): %s',hlp_tostring(scheme));
end
