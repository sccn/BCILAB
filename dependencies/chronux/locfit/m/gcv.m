function g=gcv(varargin)

fit = locfit(varargin{:});
rs = rsum(fit);

df0 = rs(1);
df1 = rs(2);
llk = rs(3);
n = size(fit.data.x,1);
gc = -2*n*llk/((n-df0)^2);
g = [llk df0 df1 gc];

return;
