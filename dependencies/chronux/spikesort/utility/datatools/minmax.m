function [extrema,inds] = minmax(X)
%MINMAX            Simultaneous overall smallest and largest components.
%   Y = MINMAX(X) returns the minimum and maximum values, of the array X
%   such that Y(1) = MIN(X) and Y(2) = MAX(X).  For N-D arrays X,
%   MINMAX(X) is equivalent to MINMAX(X(:)).
%
%   [Y,I] = MINMAX(X) also returns the linear indices of the extrema such
%   that, Y(1) == X(I(1)) and Y(2) == X(I(2)).  When X has more than one
%   extremal element, the index of the first is returned.
%
%   X must be a real, double array.  +/- Inf values behave as greater/less
%   than all other finite values, respectively.  NaN's are ignored. 
% 
%   The algorithm is significantly faster than sequentially calling
%   Matlab's min/max.  However, an already sorted input and/or a large
%   proportion of NaN values can decrease this advantage.

%%%%% Argument checking.
if (~strcmp(class(X), 'double') || ~isreal(X)), 
    error('X must be a real, double array.');  
end;

[extrema(1),extrema(2),inds(1),inds(2)] = CORE_minmax(X(:));

