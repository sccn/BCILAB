function y = sdiff(y, D, s)

%SDIFF Calculate seasonal differential of data.
%   y = SDIFF(y, D, s)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if s > 0, for i = 1:D, y = y(s+1:end)-y(1:end-s); end, end
