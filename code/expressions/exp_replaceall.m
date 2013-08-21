function [exp,changed,rules] = exp_replaceall(exp,rules)
% Apply substitution rules to some expression.
% [Out-Expression,Changed,Remaining-Rules] = exp_replaceall(In-Expression, Rule/Rule-List)
%
% Each part of the expression receives at most a single replacement, and each rule is used at most
% once rules may contain blanks and placeholders.
%
% In:
%   In-Expression : some pure expression in which substitution shall take place
%   Rule-List     : cell array of substitution rules
%
% Out:
%   Out-Expression     : In-Expression with substitutions performed
%
%   Changed            : whether a replacement has taken place
%   Remaining-Rules    : the list of rules that have not been applied
%
% Notes:
%   Impure expressions degrade into pure expressions, when substitutions are performed on them (i.e.
%   they lose their attached value). Named patterns inside the right-hand-sides of substitution
%   rules are treated as if just their symbols were given instead.
%
% See also:
%   exp_replacerepeated, exp_match
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19

if ~exp_beginfun('symbolic') return; end

[exp,changed,rules] = utl_replaceall(exp,rules);

exp_endfun;