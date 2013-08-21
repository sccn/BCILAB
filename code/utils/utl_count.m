function matches = utl_count(x,match,lev)
% Counts the number of matching subexpressions in some expression.
% Matches = utl_count(In-Expression, Match-Expression, Level)
%
% In:
%   In-Expression    : some pure expression in which subexpressions are matched
%   Match-Expression : the expression that should be matched
%   Level            : level specification, currently only Inf is allowed, meaning all levels of the expression
%
% Out:
%   Matches          : number of matching subexpressions in In-Expression
%
% See also:
%   exp_count, utl_cases
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2010-04-26


if isfield(x,{'head','parts'})
    % canonical expression: first match everything
    matches = utl_match(x,match,struct());
    % then recurse to sub-parts
    for i=1:length(x.parts)
        matches = matches + utl_count(x.parts{i},match,lev); end
elseif isfield(x,'tracking') && isfield(x.tracking,'excpression')
    % impure expression: descend
    matches = utl_count(x.tracking.expression, match, lev);
else
    % anything else...
    matches = utl_match(x,match,struct());
end
