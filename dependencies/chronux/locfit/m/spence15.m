function yhat=spence15(y)
% function for Spencer's 15-point graduation rule.
% set out following Spencer's hand-calculation method,
% which isn't the shortest computer program!

osev = ones(7,1);
n = length(y);
y = [ osev*y(1); y; osev*y(n) ];

n = length(y);
k = 3:(n-2);
a3 = y(k-1) + y(k) + y(k+1);
a2 = y(k-2) + y(k+2);
y1 = y(k)+3*(a3-a2);

n = length(y1);
k = 1:(n-3);
y2 = y1(k)+y1(k+1)+y1(k+2)+y1(k+3);

n = length(y2);
k = 1:(n-3);
y3 = y2(k)+y2(k+1)+y2(k+2)+y2(k+3);

n = length(y3);
k = 1:(n-4);
y4 = y3(k)+y3(k+1)+y3(k+2)+y3(k+3)+y3(k+4);

yhat = y4/320;
return;
