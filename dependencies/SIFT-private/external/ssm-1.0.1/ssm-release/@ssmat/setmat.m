function m = setmat(m, vec)

%@SSMAT/SETMAT Update state space matrix stationary part.
%   m = SETMAT(m, vec)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

m.mat(m.mmask)  = vec;
