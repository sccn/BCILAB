function [Y M] = ymarray(y, m, N)

%YMARRAY Create arrays of years and months.
%   [Y M] = YMARRAY(y, m, N)
%       y is the starting year.
%       m is the starting month.
%       N is the number of months.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

N1  = N - 13 + m;
if N1 > 0
    Y   = repmat(y, 1, 13-m);
    M   = m:12;
    N2  = floor(N1/12);
    Y   = [Y kron(y+1:y+N2, repmat(1, 1, 12))];
    M   = [M repmat(1:12, 1, N2)];
    N1  = mod(N1, 12);
    Y   = [Y repmat(Y(end)+1, 1, N1)];
    M   = [M 1:N1];
else
    Y   = repmat(y, 1, N);
    M   = m:m+N-1;
end
