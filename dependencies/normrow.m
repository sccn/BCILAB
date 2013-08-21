function X = normrow(X)
if nargin < 1
    error('Not enough input arguments.'); end

C = size(X,2);
if C == 1
    X = X ./ abs(X);
else
    X = sqrt(ones./(sum(X.*X,2)'))'*ones(1,C).*X;
end
