function g=lcv(varargin)
%
% likelihood cross-validation.
%
% returned vector consists of log-likelihood, trace(H), trace(H'H)
% and cross-validated deviance (-2 * cross-validated log-likelihood).
%
% Author: Catherine Loader.

fit = locfit(varargin{:},'ev','cros');
rs = rsum(fit);

df0 = rs(1);
df1 = rs(2);
llk = rs(3);

g = [llk df0 df1 -2*llk];

return;
