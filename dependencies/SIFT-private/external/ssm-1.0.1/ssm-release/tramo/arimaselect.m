function [p d q P D Q mean] = arimaselect(y, s)

%ARIMASELECT ARIMA based automatic model selection.
%   [p d q P D Q mean] = ARIMASELECT(y, s)
%   [p d q mean] = ARIMASELECT(y)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2
    % Regular ARIMA
    [d mean]    = diffdegree(y);
    y2          = y;
    if d > 0, y2 = diff(y2, d); end
    [p q]       = armadegree(y2, mean);
    P           = mean;
    fprintf(1, 'ARIMA(%d, %d, %d) ', p, d, q);
    if mean, fprintf(1, 'with mean '); end
    fprintf(1, 'selected.\n');

    arima               = ssmodel('arima', p, d, q, mean);
    [arima logL output] = estimate(y, arima, 0.1);
    fprintf(1, 'Log likelihood: %f.\nBIC: %f.\n', logL, output.BIC);
else
    % Seasonal ARIMA
    [d D mean]  = diffdegree(y, s);
    y2          = y;
    if d > 0, y2 = diff(y2, d); end
    for i = 1:D, y2 = y2(s+1:end) - y2(1:end-s); end
    [p q P Q]   = armadegree(y2, s, mean);
    fprintf(1, 'SARIMA(%d, %d, %d)*(%d, %d, %d)%d ', p, d, q, P, D, Q, s);
    if mean, fprintf(1, 'with mean '); end
    fprintf(1, 'selected.\n');

    sarima                  = ssmodel('sarima', p, d, q, P, D, Q, s, mean);
    [sarima logL output]    = estimate(y, sarima, 0.1);
    fprintf(1, 'Log likelihood: %f.\nBIC: %f.\n', logL, output.BIC);
end

