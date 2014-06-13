% [p,dp,d2p] = penPowSmooth(s,al,ep) - Smoothed power penalty
%  
% pen(s) = sqrt(s^2+ep)^al, al>0, default al=1, ep>0, default ep=1e-6
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 September 20

function [p,dp,d2p] = penPowSmooth(s,al,ep)

if nargin<2, al = 1; end                                         % default value
if nargin<3, ep = 1e-6; end                                      % default value
r   = sqrt(s.*s+ep);
p   = r.^al;                                                           % penalty
dp  = al*s.*r.^(al-2);                                        % first derivative
d2p = al*r.^(al-4).*( ep + (al-1)*s.*s);                     % second derivative
