% [p,dp,d2p] = penNegLin(s) - Linear penalty on negative part
%  
% pen(s) = max(-s,0) = -min(s,0)
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 27

function [p,dp,d2p] = penNegLin(s)

p   = max(-s,0);                                                       % penalty
dp  = -(p>0);                                                 % first derivative
d2p = zeros(size(s));                                        % second derivative