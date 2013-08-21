function [eq,dict] = exp_match(exp,form,dict)
% Check whether an expression is matched by a pattern.
% [Equals,Assignment] = exp_match(Expression, Pattern)
% 
% The pattern may contain blanks, named sub-patterns and conditional expressions. Assignment is an
% optional output which contains sub-expressions matched by named parts of the pattern.
%
% See also:
%   exp_same
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19

if ~exp_beginfun('symbolic') return; end

if ~exist('dict','var')
    dict = struct(); end

[eq,dict] = utl_match(exp,form,dict);

exp_endfun;