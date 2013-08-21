function matches = exp_count(exp,match,lev)
% Count the number of matching subexpressions in some expression (at specified levels).
% Matches = exp_count(In-Expression, Match-Expression, Level)
%
% In:
%   In-Expression    : some pure expression in which subexpressions are matched
%   Match-Expression : the expression that should be matched
%   Level            : level specification, currently only Inf is allowed, meaning all levels of the 
%                      expression
%
% Out:
%   Matches          : number of matching subexpressions in In-Expression
%
% See also:
%   exp_cases, exp_match
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-26

if ~exp_beginfun('symbolic') return; end

if ~exist('lev','var') || ~isequal(lev,Inf)
    error('currently, the only valid level is Inf, and it must be specified (for future compatibility)'); end

matches = utl_count(exp,match,lev);

exp_endfun;

