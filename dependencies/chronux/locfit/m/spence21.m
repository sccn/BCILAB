function yhat=spence21(y)
% function for Spencer's 21-point graduation rule.
% set out following Spencer's hand-calculation method,
% which isn't the shortest computer program!

oten = ones(10,1);
n = length(y);
y = [ oten*y(1); y; oten*y(n) ];

n = length(y);
k = 4:(n-3);
y1 = -y(k-3) + y(k-1) + 2*y(k) + y(k+1) -y(k+3);

n = length(y1);
k = 4:(n-3);
y2 = y1(k-3)+y1(k-2)+y1(k-1)+y1(k)+y1(k+1)+y1(k+2)+y1(k+3);

n = length(y2);
k = 3:(n-2);
y3 = y2(k-2)+y2(k-1)+y2(k)+y2(k+1)+y2(k+2);

n = length(y3);
k = 3:(n-2);
y4 = y3(k-2)+y3(k-1)+y3(k)+y3(k+1)+y3(k+2);

yhat = y4/350;
return;
