% [p,dp,d2p] = penPow(s,al) - Power penalty
%  
% pen(s) = |s|^al, al>0, default al=1
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 20

function [p,dp,d2p] = penPow(s,al)

if nargin<2, al = 1; end                                         % default value
p   = abs(s).^al;                                                      % penalty
dp  = al*sign(s).*abs(s).^(al-1);                             % first derivative
d2p = (al^2-al)*abs(s).^(al-2);                              % second derivative
