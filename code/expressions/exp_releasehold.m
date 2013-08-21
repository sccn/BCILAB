function exp = exp_releasehold(exp)
% Peel off a layer of hold expressions from some expression.
% [Out-Expression] = exp_releasehold(In-Expression)
%
% In:
%   In-Expression      : some expression
%
% Out:
%   Out-Expression     : In-Expression with the first level of Hold's peeled off
%
% Notes:
%   Descends into cell arrays.
%
% See also:
%   exp_hold
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-19

if ~exp_beginfun('symbolic') return; end

exp = utl_releasehold(exp);

exp_endfun;