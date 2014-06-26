function [alpha,b,eta,zeta,kappa] = solve_hkl_square_dual_fast(eta,data,Y,lambda,varargin);



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
tol_dx = 1e-8;
k = 1;
alpha=1;

% % optimization parameters for the dual problem
% optparam.tol = 1e-16;
% optparam.kmax = 200;
% optparam.tol_df = 1e-20;
% optparam.display= 0;
fx0 = Inf;

% solve the regular kernel learning problem
mY = mean(Y);





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

    
 
    % solve the regular kernel learning problem
    K = center_gram_matrix(double(devectorize_single( single(data.kernels * zeta ))));
    alphaK = (K + n * lambda * eye(n) ) \ ( Y - mean(Y) );
    fx =  lambda/2 * ( Y - mean(Y))'* alphaK;
    
    
    gradient_zeta = zeros(p,1);
    for i=1:p
        gradient_zeta(i) = - lambda/2 *  vectorize_quad_single(data.kernels(:,i),alphaK );
    end
    
    gradient_eta = zeros(p,1);
    for i=1:p;
        gradient_eta(i) = sum( gradient_zeta(affs{i}) .* zeta(affs{i}).^2 ) / eta(i)^2;
    end
    gradient = (1-gapeta) * gradient_eta ./ weights.^2;

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


        
    K = center_gram_matrix(double(devectorize_single( single(data.kernels * zeta ))));
    alphaK = (K + n * lambda * eye(n) ) \ ( Y - mean(Y) );
    fxalpha =  lambda/2 * ( Y - mean(Y))'* alphaK;



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

    K = center_gram_matrix(double(devectorize_single( single(data.kernels * zeta ))));
    alphaK = (K + n * lambda * eye(n) ) \ ( Y - mean(Y) );
    fxalpha =  lambda/2 * ( Y - mean(Y))'* alphaK;





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
    if norm(x-x0)<tol_dx,  break; end

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


    K = center_gram_matrix(double(devectorize_single( single(data.kernels * zeta ))));
    alpha = (K + n * lambda * eye(n) ) \ ( Y - mean(Y) );

    pred = K * alpha;
    b=mean(Y-pred );
    
kappa = zeros(p,p);
for i=1:p
    kappa(i,affs{i}) = zeta(affs{i})./eta(i);
end

 if isinf(display)
            fprintf('k=%d - f=%f -  gap=%f\n',k,fx,optcond);
        end


