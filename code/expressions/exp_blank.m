function exp = exp_blank(head)
% A blank pattern expression (for use with, e.g., exp_replace*).
% Expression = exp_blank(Head)
%
% In:
%   Head : If specified, constrains the pattern to match only expressions
%          that have a head which matches Head.
%
% Out:
%   Expression : An pattern expression, which matches any expression (with the optional requirement
%   that it has a head identical to Head).
%
% See also:
%   exp_match, exp_pattern, exp_condition, exp_rule
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-23

if exist('head','var')
    exp = struct('head',@Blank, 'parts',{{head}});
else
    exp = struct('head',@Blank, 'parts',{{}});
end