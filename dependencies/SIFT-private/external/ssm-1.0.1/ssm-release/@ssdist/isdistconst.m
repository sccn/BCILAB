function t = isdistconst(d)

%@SSDIST/ISDISTCONST True for constant non-Gaussian distributions.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

t   = ~any(d.dmask);

