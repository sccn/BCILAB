function x = x_ly(Y, M)

%X_LY Create regression variables for leap-year.
%   x = X_LY(Y, M)
%       Y is an array of years.
%       M is the corresponding array of months.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

D           = eomday(Y, M);
x           = zeros(1, length(Y));
x(D==28)    = -0.25;
x(D==29)    = 0.75;

