% [p,dp,d2p] = penAbs(s) - Absolute value penalty
%  
% pen(s) = |s|
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 20

function [p,dp,d2p] = penAbs(s)

p   = abs(s);                                                          % penalty
dp  = sign(s);                                                % first derivative
d2p = 0*s;                                                   % second derivative