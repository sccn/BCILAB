function [W,A,c,LL,Lt,gm,alpha,mu,beta,rho] = amica(x,M,m,maxiter,update_rho,mindll,iterwin,do_newton,remove_mean,As,cs)
% function [A,c,W,LL,gm,alpha,mu,beta,rho] = amica(x,[M,m,maxiter,update_rho,mindll,iterwin,do_newton,As,cs])
%
% Performa AMICA algorithm for adaptive mixutre ICA on input data x (of dimension n x N), using M models,
% with adaptive source densities unsing a mixture of m Generalized Gaussians, with fixed or adaptive shape.
%
% Stopping condition is when avg log likelihood increase is less than mindll (averaged over iterwin iterations)
% and when norm of updates is smaller than minndw.
%
% Inputs:
%
%       x           --  data, channels (n) by time points (N)

%       update_rho  --  1 = update the generalized gaussian shape parameters (default)
%                       0 = set all density shapes to 2 (Gaussian), do not update shape
%                       -rho0 = set density shapes to rho0, do not update shape

%       iterwin     --  window over which to do moving avg of log likelihood (default is 50)

%       do_newton   --  1 = use Newton update for upmixing matrices
%                       0 = use natural gradient updates

%       remove_mean --  make the data zero mean (default is 1 == yes)

%       As          --  (optional) true generating basis for plotting purposes

%       cs          --  (optional) true model centers



%
% Outputs:
%
%       W           --  learned unmixing matrices ( n x n x M )
%       A           --  learned mixing matrix with normalized columns, pinv(A) is unmixing matrix ( n x n x M )
%       c           --  learned model centers  ( n x M )
%       LL          --  history of likelihood over iterations ( 1 x min(iter,maxiter) )
%       Lt          --  log likelihood of time point for each model ( M x N )
%       gm          --  learned ica model mixture proportions ( 1 x M )
%       alpha       --  learned source density mixture proportions ( m x n x M )
%       mu          --  learned source density location parameter ( m x n x M )
%       beta        --  learned source density inverse scale parameter ( m x n x M )
%       rho         --  learned source density shape paramters ( m x n x M )
%
%
% Author: Jason Palmer (jason@sccn.ucsd.edu), Swartz Center for Computational Neuroscience, 2008.
%
% Publication: J. A. Palmer, S. Makeig, K. Kreutz-Delgado, and B. D. Rao, Newton Method for the ICA Mixture Model,
%              in Proceedings of the 33rd IEEE International Conference on Acoustics and Signal Processing (ICASSP 2008),
%              Las Vegas, NV, pp. 1805-1808, 2008.
%


if ~exist('remove_mean','var') || isempty(remove_mean)
    remove_mean = 1;
end

if nargin < 4
    maxiter = 500;
end
if nargin < 3
    m = 3;
end
if nargin < 2
    M = 1;
end
if nargin < 1
    help amica;
    return;
end

[n,N] = size(x);

lrate0 = 0.1;
lratemax = 1.0;
lnatrate = 0.1;
newt_start_iter = 5;
lratefact = 0.5;

if nargin < 6 | mindll == -1
    mindll = 1e-8;
end
if nargin < 7 | iterwin == -1
    iterwin = 50;
end

dispinterval = 50;

showLL = 1; % display the log likelihood at each iteration
retA = 1; % may cause error if W not invertible

if nargin >= 10
    plotfig = 1; % we are given the solution bases and centers
else
    plotfig = 0;  % can also set this to 1 to plot anyway (without the solution)
end

Mx = max(max(abs(x)));

rholrate = 0.5;
rhomin = 1.0;
rhomax = 10;
if nargin < 5
    update_rho = 1;
    rho0 = 1.5;
else
    if update_rho < 0
        update_rho = 0;
        rho0 = -update_rho;
    elseif update_rho == 0
        update_rho = 0;
        rho0 = 2;
    elseif update_rho == 1
        update_rho = 1;
        rho0 = 1.5;
    else
        rho0 = min(rhomax,update_rho);
        rho0 = max(rhomin,rho0);
        update_rho = 1;
    end
end

if nargin < 8
    do_newton = 1;
end

if nargin >= 10
    plothist = 1;
    nbins = 50;
else
    plothist = 0;
end

% remove the mean
if remove_mean
    for i = 1:n
        mn(i,1) = mean(x(i,:));
        x(i,:) = x(i,:) - mn(i,1);
    end
else
    mn = zeros(n,1);
end
if nargin >= 10 & max(size(cs)) > 1
    for h = 1:M
        cs(:,h) = cs(:,h) - mn;
    end
end


% initialize parameters
for h = 1:M
    W(:,:,h) = eye(n) + 0.1*randn(n,n);
    for i = 1:n
        W(i,:,h) = W(i,:,h) / norm(W(i,:,h));
    end
    c(:,h) = zeros(n,1);
end

gm = (1/M) * ones(M,1);
alpha = (1/m) * ones(m,n,M);
if m > 1
    mu = 0.1 * randn(m,n,M);
else
    mu = zeros(m,n,M);
end
beta = ones(m,n,M) + 0.1 * randn(m,n,M);

rho = rho0 * ones(m,n,M);
%rho(1,:,:) = rho0_sub * ones(1,n,M);
%rho(2:m,:,:) = rho0_sup * ones(m-1,n,M);
%rho = 2*rand(size(mu)) + 1;


% intialize variables
b = zeros(n,N,M);
g = zeros(n,N);

y = zeros(n,N,m,M);
fp = zeros(m,N);
Q = zeros(m,N);

Lt = zeros(M,N);
v = zeros(M,N);

r = zeros(n,N,m,M);

lambda = zeros(n,1);
kappa = zeros(n,1);
sigma2 = zeros(n,1);

% iteratively optimize parameters using EM
for iter = 1:maxiter
        
    % get y, Q, Lt, u (u is named z since z overwrites it later)
    for h = 1:M
        ldet(h) = log(abs(det(W(:,:,h))));
        Lt(h,:) = log(gm(h)) + ldet(h);

        if M == 1
            b = W * x;
        end
        for i = 1:n
            if M > 1
                b(i,:,h) = W(i,:,h) * x - W(i,:,h)*c(:,h);
            end
            for j = 1:m
                y(i,:,j,h) = sqrt(beta(j,i,h)) * ( b(i,:,h) - mu(j,i,h) );
                Q(j,:) = log(alpha(j,i,h)) + 0.5*log(beta(j,i,h)) + feval('logpfun',y(i,:,j,h),rho(j,i,h));
            end
            Qmax = ones(m,1)*max(Q,[],1);
            Lt(h,:) = Lt(h,:) + Qmax(1,:) + log(sum(exp(Q - Qmax),1));
            for j = 1:m
                Qj = ones(m,1)*Q(j,:);
                z(i,:,j,h) = 1 ./ sum(exp(Q-Qj),1);   % this is really u, but we use z to save space
            end
        end
    end
    Ltmax = ones(M,1) * max(Lt,[],1);
    P = sum(exp(Lt-Ltmax),1);
    LL(iter) = sum(Ltmax(1,:) + log(P)) / (n*N);
    if iter > 1
        dLL(iter) = LL(iter) - LL(iter-1);
    end
    
    if showLL
        disp(['LL(' int2str(iter) ') = ' num2str(LL(iter))]);
    end
    
    if iter > iterwin + 1
        sdll = sum(dLL(iter-iterwin+1:iter))/iterwin;
        if (sdll > 0) & (sdll < mindll)
            break
        end
        if sdll < 0
            %lrate0 = lrate0 * lratefact;
            %lratemax = lratemax * lratefact;
            %rholrate = rholrate * lratefact;
            %disp(['lrate = ' num2str(lrate)]);
        end
    end
     

    % get v, z, gm, and alpha
    for h = 1:M
        if M > 1
            Lh = ones(M,1)*Lt(h,:);
            v(h,:) = 1 ./ sum(exp(Lt-Lh),1);
            vsum(h) = sum(v(h,:));
            gm(h) = vsum(h) / N;
        
            if gm(h) == 0
                continue
            end
        end
        
        g = zeros(n,N);
        kappa = zeros(n,1);
        eta = zeros(n,1);
        sigma2 = zeros(n,1);
        
        for i = 1:n
            for j = 1:m
                if M > 1
                    z(i,:,j,h) = v(h,:) .* z(i,:,j,h);
                end
                sumz = sum(z(i,:,j,h));
                if M > 1
                    alpha(j,i,h) = sumz / vsum(h);
                else
                    alpha(j,i,h) = sumz / N;
                end

                if sumz > 0
                    z(i,:,j,h) = z(i,:,j,h) / sumz;
                else
                    continue;
                end
                fp(j,:) = feval('ffun',y(i,:,j,h),rho(j,i,h));
                zfp(j,:) = z(i,:,j,h) .* fp(j,:);
                
                g(i,:) = g(i,:) + alpha(j,i,h) * sqrt(beta(j,i,h)) * zfp(j,:);
                
                kp = beta(j,i,h) * sum(zfp(j,:).*fp(j,:));
                kappa(i) = kappa(i) + alpha(j,i,h) * kp;
                
                lambda(i) = lambda(i) + alpha(j,i,h) * (sum(z(i,:,j,h).*(fp(j,:).*y(i,:,j,h)-1).^2) + mu(j,i,h)^2 * kp);

                if rho(j,i,h) <= 2
                    if m > 1 | M > 1
                        dm = sum(zfp(j,:)./y(i,:,j,h));
                        if dm > 0
                            mu(j,i,h) = mu(j,i,h) + (1/sqrt(beta(j,i,h))) * sum(zfp(j,:)) / dm;
                        end
                    end
                    db = sum(zfp(j,:).*y(i,:,j,h));
                    if db > 0
                        beta(j,i,h) = beta(j,i,h) / db;
                    end
                else
                    if m > 1 | M > 1
                        if kp > 0
                            mu(j,i,h) = mu(j,i,h) + sqrt(beta(j,i,h)) * sum(zfp(j,:)) / kp;
                        end
                    end
                    db = ( rho(j,i,h) * sum(z(i,:,j,h).*abs(y(i,:,j,h)).^rho(j,i,h)) )^(-2 / rho(j,i,h));
                    beta(j,i,h) = beta(j,i,h) * db;
                end
                
                if update_rho
                    ytmp = abs(y(i,:,j,h)).^rho(j,i,h);
                    dr = sum(z(i,:,j,h).*log(ytmp).*ytmp);
                    dr2 = psi(1+1/rho(j,i,h)) / rho(j,i,h) - dr;
                    if ~isnan(dr2)
                        rho(j,i,h) = rho(j,i,h) + rholrate * dr2;
                    end
                    rho(j,i,h) = min(rhomax,rho(j,i,h));
                    rho(j,i,h) = max(rhomin,rho(j,i,h));    
                end                
            end
        end
        
        if M > 1
            sigma2 = b(:,:,h).^2 * v(h,:)'/vsum(h);
        else
            sigma2 = sum(b.^2,2) / N;
        end

        dW = eye(n) - g * b(:,:,h)';                                

        bflag = 0;
        for i = 1:n
            for k = 1:n            
                if i == k
                    B(i,i) = dW(i,i) / lambda(i);
                else
                    denom = kappa(i)*kappa(k)*sigma2(i)*sigma2(k) - 1;
                    if denom > 0
                        B(i,k) = (kappa(k) * sigma2(i) * dW(i,k) - dW(k,i)) / denom;
                    else
                        bflag = 1;
                    end
                end
            end
        end
        if bflag == 0 & do_newton == 1
            lrate = min(lratemax,lrate0 + (1-lrate0) * iter / newt_start_iter);
            W(:,:,h) = W(:,:,h) + lrate * B * W(:,:,h);
        else
            lrate = lnatrate;
            W(:,:,h) = W(:,:,h) + lrate * dW * W(:,:,h);
        end
    end
            
    % reparameterize
    for h = 1:M
        if gm(h) == 0
            continue
        end
        for i = 1:n
            tau = norm(W(i,:,h));
            W(i,:,h) = W(i,:,h) / tau;
            mu(:,i,h) = mu(:,i,h) / tau;
            beta(:,i,h) = beta(:,i,h) * tau^2;
        end        

        if M > 1
            % make posterior mean zero
            cnew = x * v(h,:)'/(sum(v(h,:)));
            for i = 1:n
                mu(:,i,h) = mu(:,i,h) - W(i,:,h)*(cnew-c(:,h));
            end
            c(:,h) = cnew;
        end
    end
    
    % Plot
    if plotfig & mod(iter,dispinterval) == 0
        figure(1), clf, axis(Mx*[-1 1 -1 1]);hold on;
        plot(x(1,:),x(2,:),'*');

        for h = 1:M
            A(:,:,h) = pinv(W(:,:,h));        
            for i = 1:n
                if nargin >= 10
                    plot([cs(1,h) cs(1,h)+As(1,i,h)],[cs(2,h) cs(2,h)+As(2,i,h)],'r');
                end
                cc(1) = c(1,h) + A(1,i,h)*sum(alpha(:,1,h).*mu(:,1,h));
                cc(2) = c(2,h) + A(2,i,h)*sum(alpha(:,2,h).*mu(:,2,h));
                plot([cc(1) cc(1)+A(1,i,h)],[cc(2) cc(2)+A(2,i,h)],'g');
                plot([c(1,h) c(1,h)+A(1,i,h)],[c(2,h) c(2,h)+A(2,i,h)],'y');
            end
        end
        pause(0.5);
    end
    
    if plothist & mod(iter,dispinterval) == 0
        [mval,mind] = max(Lt,[],1);
        for h = 1:M
            inds = find(mind == h);
            
            figure(2); 
            for i = 1:n
                [hst,bn] = hist(b(i,inds,h),nbins);

                dbn = abs(bn(2)-bn(1));
                hst = hst / (sum(hst)*dbn);

                ph = zeros(1,length(bn));

                subplot(M,n,n*(h-1)+i); hold off; plot(0);

                bmax = max(abs(bn(1)),abs(bn(end)));
                axis([-bmax-2*dbn bmax+2*dbn 0 max(hst)+0.1]); hold on;
                for j = 1:m
                    dh = alpha(j,i,h) * sqrt(beta(j,i,h)) * exp(-abs(sqrt(beta(j,i,h))*(bn-mu(j,i,h))).^rho(j,i,h)) / (2*gamma(1+1/rho(j,i,h)));
                    plot(bn,dh,'g');
                    ph = ph + dh;
                end
                plot(bn,hst,'r');
                plot(bn,ph);                
            end
        end
        pause(0.5);
    end

end

for h = 1:M
    if M > 1
        c(:,h) = c(:,h) + mn;  % add mean back to model centers
    end
    if retA
        for i = 1:n
            A(:,i,h) = A(:,i,h) / norm(A(:,i,h));  % normalize columns of A (mixing matrix or model basis)
        end
    end
end


function f = pfun(x,rho)
f = (1 / (2*gamma(1+1/rho))) * exp( - abs(x).^rho );

function f = logpfun(x,rho)
f = -abs(x).^rho - log(2) - gammaln(1+1/rho);

function f = ffun(x,rho)
f = rho * sign(x).*abs(x).^(rho-1);
