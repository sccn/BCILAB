function [fun grad param] = fun_cycle()

%FUN_CYCLE Create update functions for cycle component.
%   [fun grad param] = FUN_CYCLE()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

fun     = {@psi2cycle};
grad    = {@psi2cyclegrad};
param   = ssparam({'frequency'}, 'log');

function vec = psi2cycle(X)
Y       = exp(X);
sinX    = sin(Y);
cosX    = cos(Y);
vec     = [cosX; -sinX; sinX; cosX];

function [vec grad] = psi2cyclegrad(X)
Y       = exp(X);
sinX    = sin(Y);
cosX    = cos(Y);
vec     = [cosX; -sinX; sinX; cosX];
grad    = [-sinX; -cosX; cosX; -sinX];

