function [f,g] = objective(x,varargin)

% Test function for the Matlab port of LIBLBFGS
%
% Outputs: f - function value at x
%          g - gradient at x
%
% Output from the optimization should be compared with
% that from sample\sample.exe in the liblbfgs folder
%

n = numel(x);
f = 0;
g = zeros(n,1);

for i = 1:2:n-1
    t1 = 1-x(i);
    t2 = 10*(x(i+1)-x(i)*x(i));
    g(i+1) = 20*t2;
    g(i) = -2*(x(i)*g(i+1)+t1);
    f = f + t1*t1 + t2*t2;
end

end

