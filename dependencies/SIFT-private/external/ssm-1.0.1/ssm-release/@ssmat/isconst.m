function t = isconst(m)

%@SSMAT/ISCONST True for state space matrices with constant stationary part.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

t = isempty(m.mmask);

