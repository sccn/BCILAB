function x = simplex_projection(y);
n = length(y);
[ys,b]=sort(y);

% compute sequence of gradients at each point
grads(1) = sum(ys)-n*ys(1) - 1;
temp = flipud(cumsum(flipud(ys)));
% for j=2:n
%     grads(j) = temp(j) - (n-j+1) * ys(j) -1;
% end
grads = temp - 1 - ys .* ( n:-1:1)';
% find intervals where zero lies
j=min(find(grads<0))-1;
if isempty(j), j=n-1; end
% zero lies in [j,j+1];
mu = ( sum( ys(j+1:end) ) -1 ) / ( n - j );
x = max( mu - y , 0 ) + y - mu;
