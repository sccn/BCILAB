function matches = utl_cases(x,match,lev)
% Find all occurrences of some expression in another expression (at given levels).
% Matches = utl_cases(In-Expression, Match-Expression, Level)
%
% In:
%   In-Expression    : some pure expression in which subexpressions are matched
%   Match-Expression : the expression that should be matched
%   Level            : level specification, currently only Inf is allowed, meaning that all levels of the expression are checked
%
% Out:
%   Matches          : cell array of all subexpressions that are matched by Match-Expression
%
% See also:
%   exp_cases, utl_count
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2010-04-25

if all(isfield(x,{'head','parts'})) && isscalar(x)
    % canonical expression: first match everything
    if utl_match(x,match,struct())
        matches = {x};
    else
        matches = {};
    end
    % then recurse to sub-parts
    for i=1:length(x.parts)
        matches = [matches utl_cases(x.parts{i},match,lev)]; end
elseif isfield(x,'tracking') && isscalar(x) && isfield(x.tracking,'excpression')
    % impure expression: descend
    matches = utl_cases(x.tracking.expression, match, lev);
elseif utl_match(x,match,struct())
    % anything else...
    matches = {x};
else
    matches = {};
end
