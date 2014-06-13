% [p,dp,d2p] = penZero(s) - Zero penalty
%  
% pen(s) = 0
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 13

function [p,dp,d2p] = penZero(s)

p = zeros(size(s)); dp = p; d2p = p;         % penalty and first two derivatives