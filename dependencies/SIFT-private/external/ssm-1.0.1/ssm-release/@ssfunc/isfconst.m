function t = isfconst(f)

%@SSFUNC/ISFCONST True for constant state space functions.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

t   = ~any(f.fmask);

