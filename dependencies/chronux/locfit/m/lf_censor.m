function fit = lf_censor(x,y,cens,varargin)
%
% Censored local regression using normal assumption.
% Must provide x, y and cens.
% All other arguments to locfit() can be provided, with the
% exception of weights.
%
% NEED: Kaplan Meier Estimate. Iterations are fixed.
%

lfc_y = y;
unc = find(~cens);

for i = 0:3
  fit = locfit(x,lfc_y,varargin{:});
  fh = fitted(fit);

  rs = rsum(fit);
  df0 = rs(1);
  df1 = rs(2);

  rdf = sum(1-cens) - 2*df0 + df1;
  sigma = sqrt(sum( (y-fh).*(lfc_y-fh) / rdf));
  sr = (y-fh)/sigma;
  lfc_y = fh + sigma*normpdf(sr)./normcdf(-sr);
  lfc_y(unc) = y(unc);
end;

return;
