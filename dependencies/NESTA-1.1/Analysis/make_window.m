function [w,window_redundancy] = make_window(max_scale,w_type)
% [w,window_redundancy] = make_window(max_scale,w_type)
%   helper function for Peter's WFF code
%   You do NOT need to call this yourself!
%
%   w_type can be
%       'isine' for iterated sine
%       'gaussian' for gaussian
%       'trapezoid' for trapezoidal shape (not a frame)

x = (1:2^max_scale)'/2^(max_scale-1)-1;
if isequal(w_type,'isine')
    w = sin(pi/4*(1+cos(pi*x)));
    window_redundancy = 1/2;
elseif isequal(w_type,'gaussian')
    w = exp(-x.^2*8*log(2));
    window_redundancy = mean(w.^2);
elseif isequal(w_type,'trapezoid')
    w = min(1,2*(1-abs(x)));
    window_redundancy = mean(w.^2);
else
    disp('Error in make_window: unknown window type');
    disp('Options are: isine, gaussian, trapezoid');
end
