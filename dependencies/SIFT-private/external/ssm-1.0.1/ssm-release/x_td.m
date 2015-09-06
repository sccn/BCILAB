function x = x_td(Y, M, td6)

%X_TD Create regression variables for trading day.
%   x = X_TD(Y, M[, td6])
%       Y is an array of years.
%       M is the corresponding array of months.
%       Set td6 to true to use 6 variables.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, td6 = true; end

N   = length(Y);
d   = weekday(datenum(Y, M, 1));
e   = eomday(Y, M) - 28;
x   = repmat(4, 7, N);
for i = 1 : N
    x(mod((d(i):d(i)+e(i)-1)-1, 7)+1, i)    = 5;
end
if td6
    x   = x(2:7, :) - repmat(x(1, :), 6, 1);
else
    x   = sum(x(2:6, :), 1) - 2.5*sum(x([1 7], :), 1);
end

