function dsubvec = getdvec(m, mmask)

%@SSMAT/GETDVEC Get specified dynamic subvector.
%   dsubvec = GETDVEC(m, mmask)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

dsubvec = m.dvec(mmask(m.dmmask), :);

