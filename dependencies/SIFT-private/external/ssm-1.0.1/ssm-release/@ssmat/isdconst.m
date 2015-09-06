function t = isdconst(m)

%@SSMAT/ISDCONST True for state space matrices with constant dynamic part.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

t = isempty(m.dvmask);

