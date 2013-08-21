% [p,dp,d2p] = penAbsSmooth(s,ep) - Absolute value penalty
%  
% pen(s) = sqrt(s^2+ep), ep>0, default 1e-6
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 20

function [p,dp,d2p] = penAbsSmooth(s,ep)

if nargin<2, ep = 1e-6; end                                            % default
p   = sqrt(s.*s+ep);                                                   % penalty
dp  = s./p;                                                   % first derivative
d2p = ep./(p.*p.*p);                                         % second derivative