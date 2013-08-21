function z = rsum(fit)
%
% function to extract log-likelihood and degrees-of-freedom
% from locfit fit.
%
% order of returned vector:   df0 df1 llk.
%

fp = fit.fit_points;
gf = fp.kappa;

z = gf([2 3 1]);

return;
