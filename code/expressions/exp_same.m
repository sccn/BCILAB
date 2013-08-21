function res = exp_same(a,b)
% Check whether two expressions are structurally identical.
% Result = exp_same(Expression-A, Expression-B)
%
% Differs from isequal() in the following ways:
%  a) ignores the values associated with impure expressions
%  b) treats @x as equal to exp_symbol('x'); therefore, @x is a shortcut notation for exp_symbol('x').
%
% See also:
%   exp_match
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19

if ~exp_beginfun('symbolic') return; end

res = utl_same(a,b);

exp_endfun;