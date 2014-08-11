function y = arrowhead(x)
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
N = length(x);

y = zeros(N,1);
y(1) = 2*x(1)^2+sum(x.^2);
y(2:N) = x(1)^2+x(2:N).^2;
 