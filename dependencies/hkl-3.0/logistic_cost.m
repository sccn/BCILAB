function [f,grad,hess,lambda2,descent] = logistic_cost(w,x,y,lambda)
[n d] = size(x);
wTx = x * w;
b = center_logistic(wTx,y);
wTx = wTx + b;

temp = 1+exp(-wTx);
sigma = 1./temp;
 ind = find(wTx > - 50 );
f = sum( log( temp(ind) ));
inddiff=1:n;
inddiff(ind)=[];
f = f - sum( wTx(inddiff));
f = f + sum( ( 1 -y ) .* wTx ) ;
f = 1/n * f + lambda / 2 * sum(w.^2) ;
if nargout>1
	% grad = 1/n * sum(repmat( sigma - y , 1 , d ) .* x,1)' + lambda * w;
	grad = 1/n *  ( x' * ( sigma - y ) ) + lambda * w;

end
if nargout>2
    temp = sigma .* (1 - sigma );
	%hess = form_logistic_hessian(x,temp) / n +   lambda * eye(d);
	
	temp3 = repmat( 1/n * temp ,1,size(x,2)) .* x;
	hess = x' * temp3 + lambda * eye(d);
	temp2 = sum( temp3 , 1 );

	%hess = x' * ( repmat( 1/n * temp ,1,size(x,2)) .* x ) + lambda * eye(d);
	%hess = x' * sparse( diag( 1/n * temp )) * x  + lambda * eye(d);
	%hess = x' * spdiags(  1/n * temp,1,n,n ) * x  + lambda * eye(d);
	%temp2 = sum( x .* repmat( temp / n , 1, size(x,2) ) , 1 );


	hess = hess - n * (  temp2' * temp2 ) / sum(temp);
end

if nargout >3

	descent = - hess \ grad;
	lambda2 = -descent' * grad;

end
