function [x,changed,rules] = utl_replaceall(x,rules)
% Apply substitution rules to some expression.
% [Out-Expression,Changed,Remaining-Rules] = utl_replaceall(In-Expression, Rule/Rule-List)
%
% Each part of the expression receives at most a single replacement, and each rule is used at most
% once. Rules may contain blanks, named sub-expressions and conditional expressions.
%
% In:
%   In-Expression : some pure expression in which substitution shall take place
%   Rule-List     : cell array of substitution rules
%
% Out:
%   Out-Expression     : In-Expression with substitutions performed
%   Changed            : whether a replacement has taken place
%   Remaining-Rules    : the list of rules that have not been applied
%
% Notes:
%   Impure expressions degrade into pure expressions, when substitutions are performed on them (i.e.
%   they lose their attached value). Named patterns inside the right-hand-sides of substitution
%   rules are treated as if just their symbols were given instead.
%
% See also:
%   exp_replaceall, utl_replacerepeated, utl_match
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19



if isstruct(x) && isfield(x,'head') && strcmp(char(x.head),'Pattern')
	% strip off spurious pattern expressions from the expression (they are not allowed there)
    x = x.parts{1}; end

% try all rules on the full expression (which may be canonical, impure, or atomic); impurity is
% removed iff there is a match...
for r = 1:length(rules)
    [matched,dict] = utl_match(x,rules{r}.parts{1},struct());
    if matched
        % found a match; replace named sub-patterns in the substitution pattern, and return that
        if ~isempty(dict)
            x = utl_replacerepeated(rules{r}.parts{2},dict);
        else
            x = rules{r}.parts{2};
        end
        % remove the rule from the rule list
        changed = true;
        rules = rules([1:r-1 r+1:end]);
        return;
    end
end

if isfield(x,{'head','parts'})
    % the expression has substructure: try to match every subpart
    [x.head,headchanged,rules] = utl_replaceall(x.head,rules);
    partschanged = false(1,length(x.parts));
    for p=1:length(partschanged)
        [x.parts{p},partschanged(p),rules] = utl_replaceall(x.parts{p},rules); end
    changed = headchanged || any(partschanged);
else
    changed = false;
end
