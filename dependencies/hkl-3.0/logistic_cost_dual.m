function [f,grad,hess,lambda2,descent] = logistic_cost_dual(alpha,K,y,lambda)
[n] = size(K,1);

% first test if inside the domain
gamma = - n * lambda * alpha + y;
if any(gamma<1e-12) || any(gamma>1-1e-12)
	f = Inf;
	grad = NaN;
	hess = NaN;
	lambda2 = NaN;
	descent = NaN;
	return;
end

Kalpha = K * alpha;
f = lambda/2 * sum( alpha .* Kalpha ) + 1/n * sum( gamma .* log( gamma) + ( 1-gamma) .* log( 1-gamma));

if nargout>1
	grad = lambda * Kalpha - lambda * log( gamma ./ ( 1- gamma) );
end


if nargout>2
	hess = lambda * K  + lambda^2 * n * diag( 1./ ( gamma .* (1- gamma) ) );
end



if nargout >3
	% take into account the summing to zero constraint
	t1 = hess \ grad;
	t2 = hess \ ones(n,1);

	descent = - t1 + t2 * sum(t1) / sum(t2);
	lambda2 = -descent' * grad;

end
