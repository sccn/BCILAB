function y = fitted(fit)

% Fitted values from a locfit object.
%
% Input arguments:
%   fit - the locfit() fit.
%
%  Author: Catherine Loader.

z = predict(fit,'d','restyp','fit');
y = backtr(z,fit);

return;
