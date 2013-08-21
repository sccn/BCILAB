% [p,dp,d2p] = penQuad(s) - Quadratic penalty
%  
% pen(s) = s^2/2
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 13

function [p,dp,d2p] = penQuad(s)

p   = s.*s/2;                                                          % penalty
dp  = s;                                                      % first derivative
d2p = ones(size(s));                                         % second derivative