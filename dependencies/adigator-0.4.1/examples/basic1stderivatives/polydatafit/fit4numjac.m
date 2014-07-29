function p= fit4numjac(t,x, d, m)
% FIT -- Given x and d, fit() returns p
% such that norm(V*p-d) = min, where
% V = [1, x, x.?2, ... x.?(m-1)].

dim_x = size(x, 1);
if dim_x < m
  error('x must have at least m entries');
end

V = ones(dim_x, 1);

for count = 1 : (m-1)
  V = [V, x.^count];
end

p = V \ d;