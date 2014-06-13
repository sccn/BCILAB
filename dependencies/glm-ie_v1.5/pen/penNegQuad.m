% [p,dp,d2p] = penNegQuad(s) - Quadratic penalty on negative part
%  
% pen(s) = min(s,0)^2/2
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 27

function [p,dp,d2p] = penNegQuad(s)

sn  = min(s,0);                                                  % negative part 
p   = sn'*sn/2;                                                        % penalty
dp  = sn;                                                     % first derivative
d2p = double(sn<0);                                          % second derivative