function value = get(param)

%@SSPARAM/GET Get the (untransformed) parameter values.
%   value = GET(param)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

value   = [];
for i = 1 : length(param.group)-1
    X   = param.value(:, param.group(i)+1:param.group(i+1));
    switch param.transform{i}
        case '1/2 log', Y = exp(2*X);
        case 'log',     Y = exp(X);
        case 'arma1',   Y = X./realsqrt(1+X.^2);
        case 'ar2'
            Y       = X./realsqrt(1+X.^2);
            Y(:, 1) = (1-Y(:, 2)).*Y(:, 1);
        case 'ma2'
            Y       = X./realsqrt(1+X.^2);
            Y(:, 1) = (1+Y(:, 2)).*Y(:, 1);
        case 'covariance'
            p               = floor(realsqrt(2*size(X, 2)));
            Y               = zeros(size(X));
            Y(:, 1:p)       = exp(2*X(:, 1:p));
            Y(:, p+1:end)   = X(:, p+1:end)./realsqrt(1+X(:, p+1:end).^2);
            k               = p + 1;
            for i = 1:p-1, for j = i+1:p, Y(:, k) = Y(:, k).*realsqrt(Y(:, i).*Y(:, j)); k = k+1; end, end
        case 'degree of freedom', Y = 2 + exp(X);
        case 'ratio',   Y = 1./realsqrt(1+X.^2);
        case 'scale',   Y = 1 + exp(X);
        case 'l',       Y = 1.5 + 0.5*X./realsqrt(1+X.^2);
        case 'alpha',   Y = 1 + X./realsqrt(1+X.^2);
        otherwise,      Y = X; % 'identity', 'ar>=3', or 'ma>=3'
    end
    value   = [value Y];
end

