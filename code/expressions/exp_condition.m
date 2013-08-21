function res = exp_condition(pattern,condition)
% Define a conditional pattern expression, for pattern matching & substitution.
% Expression = Condition(Pattern,Condition)
%
% In:
%   Pattern : Constrains the pattern to match only expressions that match the given pattern 
%             expression; this pattern may contain named sub-expressions (e.g. using exp_pattern...)
%
%   Condition : A condition expression making use of named sub-expressions in the pattern, which 
%               needs to evaluate to true for the expression to match the conditional pattern
%
% Out:
%   An expression pattern.
%
% See also:
%   exp_pattern, exp_rule
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-23

res = struct('head',{@Condition},'parts',{{pattern,condition}});