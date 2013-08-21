% [p,dp,d2p] = penLogSmooth(s,ep) - Smoothed logarithmic penalty
%  
% pen(s) = log(s^2+ep), ep>0, default ep=1e-6
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 20

function [p,dp,d2p] = penLogSmooth(s,ep)

if nargin<2, ep = 1e-6; end                                            % default
a = s.*s+ep; p = log(a);                                               % penalty
dp  = 2*s./a;                                                 % first derivative
d2p = 2*(2*ep./a - 1)./a;                                    % second derivative