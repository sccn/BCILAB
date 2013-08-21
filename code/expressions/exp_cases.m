function matches = exp_cases(exp,match,lev)
% Find all occurrences of Match-Expression in In-Expression (at the specified levels).
% Matches = exp_cases(In-Expression, Match-Expression, Level)
%
% In:
%   In-Expression    : some pure expression in which subexpressions are matched
%   Match-Expression : the expression that should be matched
%   Level            : level specification, currently only Inf is allowed, meaning all levels of the 
%                      expression
%
% Out:
%   Matches          : cell array of all subexpressions that are matched by Match-Expression
%
% See also:
%   exp_count, exp_match, exp_same
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-25

if ~exp_beginfun('symbolic') return; end

if ~exist('lev','var') || ~isequal(lev,Inf)
    error('currently, the only valid level is Inf, and it must not be omitted (for future compatibility)'); end

matches = utl_cases(exp,match,lev);

exp_endfun;