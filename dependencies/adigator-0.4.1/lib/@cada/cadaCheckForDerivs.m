function flag = cadaCheckForDerivs(x)
% This just checks and overloaded object to see if it has any derivatives.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR

flag = 0;
for Vcount = 1:ADIGATOR.NVAROFDIFF
  if ~isempty(x.deriv(Vcount).nzlocs)
    flag = 1;
    break
  end
end