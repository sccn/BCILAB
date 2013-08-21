function res = exp_setdelayed(lhs,rhs) %#ok<STOUT>
% Delayed assignment expression; can be used in exp_block or in the global workspace.
% Rule = exp_setdelayed(Left-Hand-Side,Right-Hand-Side)
%
% The right-hand-side will be evaluated every time the left-hand-side is referenced.
%
% In:
%   Left-Hand-Side  : the target of the assignment
%   Right-Hand-Side : the expression being assigned
%
% Out:
%   Result        : the assignment expression
%
% See also:
%   exp_set
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-22

global tracking;

if ~exp_beginfun('symbolic', 'hold','all') return; end

% make an assignment in the global workspace
if ~is_symbol(lhs)
    error('For now, each lhs expression must evaluate to a symbol.'); end
tracking.stack.base.(char(get_function_symbol(lhs))) = exp_hold(rhs);

exp_endfun;