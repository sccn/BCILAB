function [exp,waschanged] = utl_replacerepeated(exp,rules)
% Apply substitution rules to some expression until it no longer changes.
% [Out-Expression,Changed] = utl_replacerepeated(In-Expression,Rule-List)
%
% rules may contain expression patterns (expressions containing blanks, named sub-expressions, conditional expressions in any mix)
%
% In:
%   In-Expression : some expression in which substitution shall take place
%   Rule-List     : cell array of substitution rules; may also be a struct, where field names are taken as symbol names
%
% Out:
%   Out-Expression : In-Expression with substitutions performed
%   Changed        : whether a replacement has taken place
%
% Notes:
%   Impure expressions degrade into pure expressions when substitutions are performed on them (i.e. they lose their attached value).
%
% See also:
%   exp_replacerepeated, utl_replaceall, utl_match
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2010-04-19

% allow for specification of rules in a struct (equivalent to an exp_rule(name,value), for every field)
if isstruct(rules) && length(rules) == 1
    dict = rules; rules = {};
    for fn = fieldnames(dict)'
        rules{end+1} = exp_rule(exp_symbol(fn{1}),dict.(fn{1})); end %#ok<AGROW>
end

waschanged = false;
while 1
    % perform replacements until the expression no longer changes...
    [exp,changed] = utl_replaceall(exp,rules);
    if changed
        waschanged = true;
    else
        break;        
    end
end
