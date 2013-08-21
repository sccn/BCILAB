function res = exp_rule(match,repl)
% Create a substitution rule, which to replace a Pattern in an expression.
% Rule = exp_rule(Pattern,Replacement)
%
% In:
%   Pattern     : An expression; may contain blanks (see exp_blank()), named sub-expressions (see
%                 exp_pattern()) or generally any structure (see exp_symbol()).
%   Replacement : a replacement expression; may contain symbols and/or named patterns that reference 
%                 parts inside Pattern.
%
% Out:
%   Rule        : a replacement rule
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-23

res = struct('head',{@Rule},'parts',{{match,repl}});