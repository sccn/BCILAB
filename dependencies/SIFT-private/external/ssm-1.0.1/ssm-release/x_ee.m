function x = x_ee(Y, M, d)

%X_EE Create regression variables for Easter effect.
%   x = X_EE(Y, M, d)
%       Y is an array of years.
%       M is the corresponding array of months.
%       d is the number of days before Easter.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

AVG_RATIO3  = 0.35332;%0.232;
AVG_RATIO4  = 0.64668;%0.768;

x           = zeros(1, length(Y));

% March
M3          = (M == 3);
[EM3 ED3]   = easter(Y(M3));

M33         = M3;
M33(M33)    = (EM3 == 3);
x(M33)      = 1 - AVG_RATIO3;

M34         = M3;
M34(M34)    = (EM3 == 4);
RATIO3      = (d + 1 - ED3(EM3 == 4))/d;
RATIO3(RATIO3<0)    = 0;
x(M34)      = RATIO3 - AVG_RATIO3;

% April
M4          = (M == 4);
[EM4 ED4]   = easter(Y(M4));

M43         = M4;
M43(M43)    = (EM4 == 3);
x(M43)      = -AVG_RATIO4;

M44         = M4;
M44(M44)    = (EM4 == 4);
RATIO4      = (ED4(EM4 == 4) - 1)/d;
RATIO4(RATIO4>1)    = 1;
x(M44)      = RATIO4 - AVG_RATIO4;

