function b = center_logistic(u,y);
s = mean(y);
b = log( s / (1 - s ) ) - mean(u);
for i=1:20
	sigma = 1./ ( 1 + exp( - u - b ) );
	phider = sum(sigma - y) ;
	phider2 = sum(sigma .* ( 1- sigma));
	bnew = b - phider / phider2;
	if norm(bnew-b) < 1e-16, break; end
	b=bnew;

end

