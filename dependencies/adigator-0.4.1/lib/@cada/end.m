function out = end(x,dim,ndim)
% CADA overloaded END
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
if ndim == 1
  out = length(x);
else
  out = size(x,dim);
end