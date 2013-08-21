function res = exp_set(lhs,rhs) %#ok<STOUT>
% An assignment expression; can be used inside exp_block or against the global workspace.
% Rule = exp_set(Left-Hand-Side,Right-Hand-Side)
%
% In:
%   Left-Hand-Side  : the target of the assignment
%   Right-Hand-Side : the expression being assigned
%
% Out:
%   Result        : the assignment expression
%
% See also:
%   exp_setdelayed
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-22

global tracking;

if ~exp_beginfun('symbolic', 'hold','first') return; end

% make an assignment in the global workspace
if ischar(lhs)
    tracking.stack.base.(lhs) = rhs;
else
    if ~is_symbol(lhs)
        error('For now, each lhs expression must evaluate to a symbol.'); end
    tracking.stack.base.(char(get_function_symbol(lhs))) = rhs;
end

exp_endfun;
