function g=aic(varargin)

pen = 2.0;

fit = locfit(varargin{:});
rs = rsum(fit);

df0 = rs(1);
df1 = rs(2);
llk = rs(3);
ai = -2*llk + pen * df0;
g = [llk df0 df1 ai];

return;
