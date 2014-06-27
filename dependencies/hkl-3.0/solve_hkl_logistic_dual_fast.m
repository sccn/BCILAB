function [alpha,b,eta,zeta,kappa] = solve_hkl_logistic_dual_fast(eta,data,Y,lambda,varargin);



% solve the milasso for a given regularization parameter
n = size(data.X,1);
p =length(data.affinity);
pX = size(data.X,2);
affs = data.affinity;
weights = data.weights;

% optional parameters
mingap = 1e-3;
display = 0;
gapeta = 1e-3;

% READ OPTIONAL PARAMETERS
args = varargin;
nargs = length(args);
for i=1:2:nargs
	switch args{i},
		case 'mingap',        mingap = args{i+1};
		case 'gapeta',        gapeta = args{i+1};
		case 'display',        display = args{i+1};
	end
end


if isempty(eta)
	eta = 1/p ./ weights.^2;
end
x = eta .* weights.^2;
x = max(x,0);
x = x / sum(x);

if display==1
	display = Inf;
end
beta_ARMIJO = .25 ;
sigma_ARMIJO = .1 ;
kmax_ARMIJO = 20;
kmax = 200;
tol_dx = 1e-12;
k = 1;
alpha=1;

% % optimization parameters for the dual problem
% optparam.tol = 1e-16;
% optparam.kmax = 200;
% optparam.tol_df = 1e-20;
% optparam.display= 0;
fx0 = Inf;

% solve the regular kernel learning problem
X = data.X;	% no need to center here for logistic
mY = mean(Y);
G = zeros(size(data.X,2),p);
for i=1:p
	G(data.groups{i},i) = 1;
end


% optimization parameter for logistic regression
optparam.tol = 1e-16;
optparam.kmax = 50;
optparam.tol_df = 1e-14;
optparam.display= 0;


%nu =  1e-6; % added ridge on kernels
while k<=kmax;
	x = max(x,0);
	x = x / sum(x);


	% compute value of function
	eta = ( x*(1-gapeta) + gapeta/p )./ weights.^2 ;
	zeta = zeros(p,1);
	for i=1:p
		zeta(affs{i}) = zeta(affs{i}) + 1./eta(i);
	end
	zeta = 1./zeta;
		reweight = G * sqrt(zeta);

	% solve the regular kernel learning problem in the primal
	if k==1
		% first iteration -> initialize with square loss
		% otherwise, take previous iteration weights
		temp = ( reweight .* ( X' * (Y - mY ) ) );
		wt = (  (reweight*reweight') .* ( X' * X ) + n * lambda * eye(pX) ) \   temp ;
	end

	Xr = repmat(reweight',n,1) .* X;
	[wt,exit_parameters] = minimize_newton(wt,@logistic_cost,optparam, Xr,Y,lambda);
	fx = logistic_cost(wt,Xr,Y,lambda);
 

	gradient_zeta = zeros(p,1);
	for i=1:p
		gradient_zeta(i) = - lambda/2 *  norm(wt(data.groups{i}))^2 / zeta(i);
	end

	gradient_eta = zeros(p,1);
	for i=1:p;
		gradient_eta(i) = sum( gradient_zeta(affs{i}) .* zeta(affs{i}).^2 ) / eta(i)^2;
	end
	gradient = (1-gapeta) * gradient_eta ./ weights.^2;


 alpha = 1;
	% check optimality condition
	optcond = max(x .* ( gradient - min(gradient)));
	if max(x .* ( gradient - min(gradient))) < mingap, break; end

	fx0=fx;

	% projected gradient
	xbar = x -   gradient;
	d = - gradient;
	x0=x;
	ialpha = 1;
	while ialpha<kmax_ARMIJO
		xalpha = simplex_projection(x + d * alpha);
 		% fxalpha = feval(fun,xalpha,varargin{:});
		eta = ( xalpha*(1-gapeta) + gapeta/p )./ weights.^2 ;
		zeta = zeros(p,1);
		for i=1:p
			zeta(affs{i}) = zeta(affs{i}) + 1./eta(i);
		end
		zeta = 1./zeta;



		reweight = G * sqrt(zeta);
		Xr = repmat(reweight',n,1) .* X;
		[wt0,exit_parameters] = minimize_newton(wt,@logistic_cost,optparam, Xr,Y,lambda);
		fxalpha = logistic_cost(wt0,Xr,Y,lambda);

 		if fx - fxalpha >= - sigma_ARMIJO * gradient' * ( xalpha - x );
 			% start to go up
			while ( fx - fxalpha >= - sigma_ARMIJO * gradient' * ( xalpha - x ) ) & ialpha<kmax_ARMIJO;
				alpha = alpha / beta_ARMIJO;
 				xalpha = simplex_projection(x + d * alpha);
				%fxalpha = feval(fun,xalpha,varargin{:});
				eta = ( xalpha*(1-gapeta) + gapeta/p )./ weights.^2 ;
				zeta = zeros(p,1);
				for i=1:p
					zeta(affs{i}) = zeta(affs{i}) + 1./eta(i);
				end
				zeta = 1./zeta;

				reweight = G * sqrt(zeta);
				Xr = repmat(reweight',n,1) .* X;
				[wt0,exit_parameters] = minimize_newton(wt,@logistic_cost,optparam, Xr,Y,lambda);
				fxalpha = logistic_cost(wt0,Xr,Y,lambda);





				% fprintf('+');
				ialpha = ialpha+1;
			end
			alpha = alpha * beta_ARMIJO;
			xalpha = simplex_projection(x + d * alpha);
			x = xalpha;
			break;
		else
			alpha = alpha * beta_ARMIJO;
 		end
		ialpha = ialpha+1;
	end
 	if ialpha==kmax_ARMIJO, alpha=1; end



	if ~isinf(display)
		if display==1,
			fprintf('k=%d - f=%f - armijo=%d - dx=%e - gap=%f\n',k,fx,ialpha,norm(x-x0),optcond);
		end
		if display>1,
			if mod(k,display)==1
				fprintf('k=%d - f=%f - armijo=%d - dx=%e - gap=%f\n',k,fx,ialpha,norm(x-x0),optcond);

			end
		end
	end
	if norm(x-x0)<tol_dx,    break; end

	if sum(x) > 1+1e-1, keyboard; end

	k=k+1;
end


% if k==kmax+1
%     keyboard;
% end


eta = ( x*(1-gapeta) + gapeta/p )./ weights.^2 ;
zeta = zeros(p,1);
for i=1:p
	zeta(affs{i}) = zeta(affs{i}) + 1./eta(i);
end
zeta = 1./zeta;
reweight = G * sqrt(zeta);
Xr = repmat(reweight',n,1) .* X;
[wt,exit_parameters] = minimize_newton(wt,@logistic_cost,optparam, Xr,Y,lambda);

pred = data.X * ( reweight.* wt);
b= center_logistic(pred,Y);

alpha = 1/(n*lambda) * ( Y - 1./( 1 + exp( - pred - b ) ) );

kappa = zeros(p,p);
for i=1:p
	kappa(i,affs{i}) = zeta(affs{i})./eta(i);
end

if isinf(display)
	fprintf('k=%d - f=%f -  gap=%f\n',k,fx,optcond);
end
% if optcond> 1,
% keyboard
% end

