function z = dot(x,y,varargin)
% CADA overloaded DOT - just calls z = sum(x.*y) - lazy.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
if nargin == 2
  z = sum(x.*y);
else
  z = sum(x.*y,varargin{1});
end