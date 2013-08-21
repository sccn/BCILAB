function y = residuals(fit,type)

% Residuals (or a few other things) from a locfit() fit.
%
% Input arguments:
%   fit - the locfit() fit.
%   type (optional) type of residuals. Valid types are
%     'dev'    (deviance, the default)
%     'd2'     (deviance squared)
%     'pearson'(Pearson)
%     'raw'    (observed - fitted)
%     'ldot'   (derivative of log-likelihood)
%     'lddot'  (second derivative)
%     'fit'    (fitted values - no transformation)
%     'mean'   (fitted values - with back transformation)
%
%  Author: Catherine Loader.

if (nargin<2) type = 'dev'; end;
  
y = predict(fit,'d','restyp',type);

return;
