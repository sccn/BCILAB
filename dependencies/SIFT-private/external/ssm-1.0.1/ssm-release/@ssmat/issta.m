function t = issta(m)

%@SSMAT/ISSTA True for stationary state space matrices.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

t = isempty(m.dmmask);

