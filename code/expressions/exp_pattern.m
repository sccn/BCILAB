function res = exp_pattern(symb,expr)
% Define a named pattern expression, optionally with a specific structure.
% Expression = Pattern(Symbol,Pattern)
%
% In:
%   Symbol  : A symbol, expressed as @name or exp_symbol('name')
%
%   Pattern : Constrains the pattern to match only expressions that match the given pattern
%             expression. Defaults to exp_blank(), which matches any expression.
%
% Out:
%   An expression pattern, which matches any or the given pattern expression; can be used, among 
%   others, in exp_replace* to refer to specific parts of an expression.
%
% See also:
%   exp_declare_patterns
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-28

if ~is_symbol(symb)
    error('Only symbols may be specified as pattern names.'); end
if ~exist('expr','var')
    expr = exp_blank(); end

res = struct('head',{@Pattern},'parts',{{symb,expr}});