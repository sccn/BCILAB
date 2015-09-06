function param = set(param, value)

%@SSPARAM/SET Set the (untransformed) parameter values.
%   param = SET(param, value)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

tol         = 100*eps;

param.value = [];
for i = 1 : length(param.group)-1
    Y   = value(:, param.group(i)+1:param.group(i+1));
    switch param.transform{i}
        case '1/2 log', X = reallog(Y)/2;
        case 'log',     X = reallog(Y);
        case 'arma1'
            D   = 1 - Y.^2;
            if any(D < 0 & abs(D) < tol), X = Y*Inf; warning('ssm:ssparam_set:RoundingError', 'Rounding error (?) detected and corrected.');
            else X = Y./realsqrt(D); end
        case 'ar2'
            X       = zeros(size(Y));
            D       = 1 - Y(:, 2).^2;
            if any(D < 0 & abs(D) < tol), X(:, 2) = Y(:, 2)*Inf; warning('ssm:ssparam_set:RoundingError', 'Rounding error (?) detected and corrected.');
            else X(:, 2) = Y(:, 2)./realsqrt(D); end
            X(:, 1) = Y(:, 1)./(1-Y(:, 2));
            D       = 1 - X(:, 1).^2;
            if any(D < 0 & abs(D) < tol), X(:, 1) = X(:, 1)*Inf; warning('ssm:ssparam_set:RoundingError', 'Rounding error (?) detected and corrected.');
            else X(:, 1) = X(:, 1)./realsqrt(D); end
        case 'ma2'
            X       = zeros(size(Y));
            D       = 1 - Y(:, 2).^2;
            if any(D < 0 & abs(D) < tol), X(:, 2) = Y(:, 2)*Inf; warning('ssm:ssparam_set:RoundingError', 'Rounding error (?) detected and corrected.');
            else X(:, 2) = Y(:, 2)./realsqrt(D); end
            X(:, 1) = Y(:, 1)./(1+Y(:, 2));
            D       = 1 - X(:, 1).^2;
            if any(D < 0 & abs(D) < tol), X(:, 1) = X(:, 1)*Inf; warning('ssm:ssparam_set:RoundingError', 'Rounding error (?) detected and corrected.');
            else X(:, 1) = X(:, 1)./realsqrt(D); end
        case 'covariance'
            p               = floor(realsqrt(2*size(Y, 2)));
            X               = zeros(size(Y));
            X(:, 1:p)       = reallog(Y(:, 1:p))/2;
            k               = p + 1;
            for i = 1:p-1, for j = i+1:p, X(:, k) = Y(:, k)./realsqrt(Y(:, i).*Y(:, j)); k = k+1; end, end
            X(:, p+1:end)   = X(:, p+1:end)./realsqrt(1-X(:, p+1:end).^2);
        case 'degree of freedom', X = reallog(Y-2);
        case 'ratio',   X = realsqrt(1./Y.^2-1);
        case 'scale',   X = reallog(Y-1);
        case 'l',       X = 2*(Y-1.5)./realsqrt(1-4*(Y-1.5).^2);
        case 'alpha',   X = (Y-1)./realsqrt(1-(Y-1).^2);
        otherwise,      X = Y; % 'identity', 'ar>=3', or 'ma>=3'
    end
    param.value = [param.value X];
end

