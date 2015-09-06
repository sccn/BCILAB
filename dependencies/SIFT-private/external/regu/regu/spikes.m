function [A,b,x] = spikes(n,t_max)
%SPIKES Test problem with a "spiky" solution.
%
% [A,b,x] = spikes(n,t_max)
%
% Artificially generated discrete ill-posed problem.
%
% The solution x consists of a unit step at t = .5, and a pulse train
% of spikes of decrasing magnitude at t = .5, 1.5, 2.5, ...
%
% The parameter t_max is optional; its default value is 5.
% It controls the length of the pulse train.

% Per Christian Hansen, IMM, 04/21/97.

% Initialization.
if (nargin == 1), t_max = 5; end
del = t_max/n;

% Compute the matrix A.
[t,sigma] = meshgrid(del:del:t_max,del:del:t_max);
A = sigma./(2*sqrt(pi*t.^3)).*exp(-(sigma.^2)./(4*t));

% Compute the right-hand side b and the solution x.
if (nargout > 1)
  heights = 2*ones(t_max,1); heights(1) = 25;
  heights(2) = 9; heights(3) = 5; heights(4) = 4; heights(5) = 3;
  x = zeros(n,1); n_h = 1;
  peak = 0.5/t_max; peak_dist = 1/t_max;
  if (peak < 1)
    n_peak = round(peak*n); x(n_peak) = heights(n_h);
    x(n_peak+1:n) = ones(n-n_peak,1);
    peak = peak + peak_dist; n_h = n_h + 1;
  end
  while (peak < 1)
    x(round(peak*n)) = heights(n_h);
    peak = peak + peak_dist; n_h = n_h + 1;
  end
  b = A*x;
end