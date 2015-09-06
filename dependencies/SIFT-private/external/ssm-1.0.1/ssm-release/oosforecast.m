function [yf err SS estmodels] = oosforecast(y, model, n1, h, varargin)

%OOSFORECAST Out-of-sample forecast and diagnostics.
%   [yf err SS] = OOSFORECAST(y, model, n1, h)
%       n1 is the number of time points to exclude at the end.
%       h is the number of steps ahead to forecast, can be an array.
%       yf is the forecast obtained.
%       err is the forecast error.
%       SS is the forecast error cumulative sum of squares.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

minh    = min(h);
p       = size(y, 1);
n       = size(y, 2);
yf      = repmat(NaN, [p n1-minh+1 length(h)]);
param0  = model.param;
Zdyn    = ~issta(model.Z);
Zmat    = getmat(model.Z);
if ~Zdyn, Z = Zmat; end
if nargout > 3, estmodels = cell(1, n1-minh+1); end
for t = n-n1 : n-minh
    model           = estimate(y(:, 1:t), model, param0, [], varargin{:});
    afore           = kalman([y(:, 1:t) repmat(NaN, p, n-t)], model, varargin{:});
    for i = 1 : length(h)
        t1  = t+h(i);
        if t1 <= n
            if Zdyn, Z = Zmat{t1}; end
            yf(:, t1-n+n1-minh+1, i)    = Z*afore(:, t1);
        end
    end
    if nargout > 3, estmodels{t-n+n1+1} = model; end
end

err             = yf - repmat(y(:, n-n1+minh:n), [1 1 length(h)]);
SS              = err;
SS(isnan(SS))   = 0;
SS              = cumsum(SS.^2, 2);
SS(isnan(err))  = NaN;

if p == 1
    yf      = permute(yf, [3 2 1]);
    err     = permute(err, [3 2 1]);
    SS      = permute(SS, [3 2 1]);
end
