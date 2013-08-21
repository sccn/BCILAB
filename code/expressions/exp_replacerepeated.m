function [exp,waschanged] = exp_replacerepeated(exp,rules)
% Apply substitution rules to some expression until it no longer changes.
% [Out-Expression,Changed] = exp_replacerepeated(In-Expression,Rule-List)
%
% Rules may contain blanks and placeholders.
%
% In:
%   In-Expression : some expression in which substitution shall take place
%   Rule-List     : cell array of substutition rules; may also be a struct, where field names are 
%                   taken as symbol names
%
% Out:
%   Out-Expression     : In-Expression with substitutions performed
%   Changed            : whether a replacement has taken place
%
% Notes:
%   Impure expressions degrade into pure expressions when substitutions are performed on them (i.e.
%   they lose their value).
%
% See also:
%   exp_replaceall, exp_match
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19

if ~exp_beginfun('symbolic') return; end

[exp,waschanged] = utl_replacerepeated(exp,rules);

exp_endfun;