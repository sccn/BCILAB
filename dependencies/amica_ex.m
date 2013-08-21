function [A,c,LL,Lt,gm,alpha,mu,beta,rho] = amica_ex(x,M,m,maxiter,update_rho,mindll,iterwin,do_newton,remove_mean,As,cs)
% function [A,c,LL,gm,alpha,mu,beta,rho] = amica_ex(x,[M,m,maxiter,update_rho,mindll,iterwin,do_newton,As,cs])
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
%       M           --  number of ICA mixture models (default is 1)
%       m           --  number of source density mixtures (default is 3)
%       maxiter     --  maximum number of iterations to perform (default is 500)
%       update_rho  --  1 = update the generalized gaussian shape parameters (default)
%                       0 = set all density shapes to 2 (Gaussian), do not update shape
%                       -rho0 = set density shapes to rho0, do not update shape
%       mindll      --  stop when avg log likelihood changing less than this (default is 1e-8)
%       iterwin     --  window over which to do moving avg of log likelihood (default is 50)
%       do_newton   --  1 = use Newton update for upmixing matrices
%                       0 = use natural gradient updates
%       remove_mean --  make the data zero mean (default is 1 == yes)
%       As          --  (optional) true generating basis for plotting purposes
%       cs          --  (optional) true model centers
%
% Outputs:
%
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

if nargin < 4
    maxiter = 100;
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
if nargin >= 10 && nargin < 11
    cs = zeros(n,M);
end

lrate0 = 0.1;
lratemax = 1.0;
lnatrate = 0.1;
lnatratemax = 0.1;
newt_start_iter = 25;
lratefact = 0.5;

maxdec = 3;
numdec = 0;
maxdec2 = 20;
numdec2 = 0;

lrate = lrate0;
remove_mean = 1;

if nargin < 6 | mindll == -1
    mindll = 1e-8;
end
if nargin < 7 | iterwin == -1
    iterwin = 1;
end

if M > 1
    dispinterval = 10;
else
    dispinterval = 1;
end

if nargin >= 10
    showLL = 1; % display the log likelihood at each iteration
else
    showLL = 1;
end
if nargin >= 10
    plotfig = 1; % we are given the solution bases and centers
else
    plotfig = 0;  % can also set this to 1 to plot anyway (without the solution)
end

Mx = max(max(abs(x)));

rholrate = 0.1;
rhomin = 1.0;
rhomax = 2;
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
if nargin < 9
    remove_mean = 1;
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
    A(:,:,h) = eye(n) + 0.1*randn(n,n);
    for i = 1:n
        A(:,i,h) = A(:,i,h) / norm(A(:,i,h));
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
v = ones(M,N);
z = ones(n,N,m,M)/N;

r = zeros(n,N,m,M);

lambda = zeros(n,1);
kappa = zeros(n,1);
sigma2 = zeros(n,1);

% iteratively optimize parameters using EM
for iter = 1:maxiter
        
    % get y, Q, Lt, u (u is named z since z overwrites it later)
    for h = 1:M
        ldet(h) = -log(abs(det(A(:,:,h))));
        Lt(h,:) = log(gm(h)) + ldet(h);

        if M == 1
            b = pinv(A) * x;
        end
        for i = 1:n
            if M > 1
                Wh = pinv(A(:,:,h));
                b(i,:,h) = Wh(i,:) * x - Wh(i,:)*c(:,h);
            end
            for j = 1:m
                y(i,:,j,h) = sqrt(beta(j,i,h)) * ( b(i,:,h) - mu(j,i,h) );
                Q(j,:) = log(alpha(j,i,h)) + 0.5*log(beta(j,i,h)) + feval('logpfun',y(i,:,j,h),rho(j,i,h));
            end
            if m > 1
                Qmax = ones(m,1)*max(Q,[],1);
                Lt(h,:) = Lt(h,:) + Qmax(1,:) + log(sum(exp(Q - Qmax),1));
                for j = 1:m
                    Qj = ones(m,1)*Q(j,:);
                    z(i,:,j,h) = 1 ./ sum(exp(Q-Qj),1);   % this is really u, but we use z to save space
                end
            else
                Lt(h,:) = Lt(h,:) + Q(1,:);
            end
        end
    end
    if M > 1
        Ltmax = ones(M,1) * max(Lt,[],1);
        P = sum(exp(Lt-Ltmax),1);
        LL(iter) = sum(Ltmax(1,:) + log(P)) / (n*N);
    else
        LL(iter) = sum(Lt) / (n*N);
    end
        
    if iter > 1
        dLL(iter) = LL(iter) - LL(iter-1);
    end
    
    if showLL
        disp(['iter '  int2str(iter) ' lrate = ' num2str(lrate) ' LL = ' num2str(LL(iter))]);
    end
    
    if iter > iterwin + 1
        sdll = sum(dLL(iter-iterwin:iter))/iterwin;
        if (sdll > 0) && (sdll < mindll)
            break
        end
        if sdll < 0
            disp('Likelihood decreasing!');
            lrate = lrate * lratefact;
            numdec = numdec + 1;
            if numdec > maxdec
                disp(['Reducing lrate, lrate = ' num2str(lrate)]);
                lratemax = lratemax * lratefact;
                lnatratemax = lnatratemax * lratefact;
                rholrate = rholrate * lratefact;
                numdec2 = numdec2 + 1;
                if numdec2 > maxdec2
                    break;
                end
                numdec = 0;
            end
        else
            if iter > newt_start_iter && do_newton
                lrate = min(lratemax,lrate + min(0.1,lrate));
            else
                lrate = min(lnatratemax,lrate + min(0.1,lrate));
            end
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
        lambda = zeros(n,1);
        eta = zeros(n,1);
        sigma2 = zeros(n,1);
        
        for i = 1:n
            for j = 1:m
                if M > 1
                    if m > 1
                        z(i,:,j,h) = v(h,:) .* z(i,:,j,h);
                        sumz = sum(z(i,:,j,h));
                        alpha(j,i,h) = sumz / vsum(h);
                    else
                        z(i,:,j,h) = v(h,:);
                        sumz = sum(z(i,:,j,h));
                    end
                else
                    if m > 1
                        sumz = sum(z(i,:,j,h));
                        alpha(j,i,h) = sumz / N;
                    else
                        sumz = N;
                    end
                end

                if sumz > 0
                    if M > 1 | m > 1
                        z(i,:,j,h) = z(i,:,j,h) / sumz;
                    end
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
                    if rho(j,i,h) > 2
                        dr2 = psi(1+1/rho(j,i,h)) / rho(j,i,h) - dr;
                        if ~isnan(dr2)
                            rho(j,i,h) = rho(j,i,h) + 0.5 * dr2;
                        end
                    else
                        dr2 = 1 - rho(j,i,h) * dr / psi(1+1/rho(j,i,h));
                        if ~isnan(dr2)
                            rho(j,i,h) = rho(j,i,h) + rholrate * dr2;
                        end
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

        dA = eye(n) - g * b(:,:,h)';

        bflag = 0;
        for i = 1:n
            for k = 1:n            
                if i == k
                    B(i,i) = dA(i,i) / (-0*dA(i,i) + lambda(i));
                else
                    denom = kappa(i)*kappa(k)*sigma2(i)*sigma2(k) - 1;
                    if denom > 0
                        B(i,k) = (-kappa(k) * sigma2(i) * dA(i,k) + dA(k,i)) / denom;
                    else
                        bflag = 1;
                    end
                end
            end
        end
        if bflag == 0 && do_newton == 1 && iter > newt_start_iter
            A(:,:,h) = A(:,:,h) + lrate * A(:,:,h) * B;
        else
            A(:,:,h) = A(:,:,h) - lnatrate * A(:,:,h) * dA;
        end
    end
            
    % reparameterize
    for h = 1:M
        if gm(h) == 0
            continue
        end
        for i = 1:n
            tau = norm(A(:,i,h));
            A(:,i,h) = A(:,i,h) / tau;
            mu(:,i,h) = mu(:,i,h) * tau;
            beta(:,i,h) = beta(:,i,h) / tau^2;
        end

        if M > 1
            % make posterior mean zero
            cnew = x * v(h,:)'/(sum(v(h,:)));
            for i = 1:n
                Wh = pinv(A(:,:,h));
                mu(:,i,h) = mu(:,i,h) - Wh(i,:)*(cnew-c(:,h));
            end
            c(:,h) = cnew;
        end
    end
    
    % Plot
    if plotfig & mod(iter,dispinterval) == 0
        if M < 2
            figure(1), clf;
            subplot(2,2,1); set(gca,'fontsize',14); hold on; axis(Mx*[-1 1 -1 1]);
            xlabel('\itx_{\rm1}'); ylabel('\itx_{\rm2}');
            plot(x(1,:),x(2,:),'g.');
            
            for h = 1:M
                for i = 1:n
                    if nargin >= 10
                        plot([cs(1,h) cs(1,h)+As(1,i,h)],[cs(2,h) cs(2,h)+As(2,i,h)],'r','LineWidth',2);
                    end
                    %cc(1) = c(1,h) + A(1,i,h)*sum(alpha(:,1,h).*mu(:,1,h));
                    %cc(2) = c(2,h) + A(2,i,h)*sum(alpha(:,2,h).*mu(:,2,h));
                    %plot([cc(1) cc(1)+A(1,i,h)],[cc(2) cc(2)+A(2,i,h)],'b','LineWidth',2);
                    plot([c(1,h) c(1,h)+A(1,i,h)],[c(2,h) c(2,h)+A(2,i,h)],'b','LineWidth',2);
                end
            end
            %pause(0.5);

            [mval,mind] = max(Lt,[],1);
            for h = 1:1
                inds = find(mind == h);
                
                for i = 1:2
                    [hst,bn] = hist(b(i,inds,h),nbins);
                    
                    dbn = abs(bn(2)-bn(1));
                    hst = hst / (sum(hst)*dbn);
                    
                    ph = zeros(1,length(bn));
                    
                    subplot(2,2,1+i); set(gca,'fontsize',14); hold on;
                    xlabel(['\its_{\rm' int2str(i) '}']); ylabel('norm. histogram');
                    
                    bmax = max(abs(bn(1)),abs(bn(end)));
                    axis([-bmax-2*dbn bmax+2*dbn 0 max(hst)+0.1]); hold on;
                    for j = 1:m
                        dh = alpha(j,i,h) * sqrt(beta(j,i,h)) * exp(-abs(sqrt(beta(j,i,h))*(bn-mu(j,i,h))).^rho(j,i,h)) / (2*gamma(1+1/rho(j,i,h)));
                        plot(bn,dh,'g');
                        ph = ph + dh;
                    end
                    plot(bn,hst,'r');
                    plot(bn,ph);
                    axis([-Mx Mx 0 max(max(ph,1))])
                end
            end
            subplot(2,2,4); set(gca,'fontsize',14); plot(LL); xlim([1 maxiter]);
            xlabel('iteration'); ylabel('Log Likelihood');
        else
            figure(1), clf;
            subplot(1,2,1); set(gca,'fontsize',14); hold on; axis(Mx*[-1 1 -1 1]);
            xlabel('\itx_{\rm1}'); ylabel('\itx_{\rm2}');
            plot(x(1,:),x(2,:),'g.');
            
            for h = 1:M
                for i = 1:n
                    if nargin >= 10
                        plot([cs(1,h) cs(1,h)+As(1,i,h)],[cs(2,h) cs(2,h)+As(2,i,h)],'r','LineWidth',2);
                    end
                    %cc(1) = c(1,h) + A(1,i,h)*sum(alpha(:,1,h).*mu(:,1,h));
                    %cc(2) = c(2,h) + A(2,i,h)*sum(alpha(:,2,h).*mu(:,2,h));
                    %plot([cc(1) cc(1)+A(1,i,h)],[cc(2) cc(2)+A(2,i,h)],'b','LineWidth',2);
                    plot([c(1,h) c(1,h)+A(1,i,h)],[c(2,h) c(2,h)+A(2,i,h)],'b','LineWidth',2);
                end
            end
            subplot(1,2,2); set(gca,'fontsize',14); plot(LL); xlim([1 maxiter]);
            xlabel('iteration'); ylabel('Log Likelihood');
            
        end
    end
end

for h = 1:M
    if M > 1
        c(:,h) = c(:,h) + mn;  % add mean back to model centers
    end
end
disp('Done!');

function f = pfun(x,rho)
f = (1 / (2*GAMMA(1+1/rho))) * exp( - abs(x).^rho );

function f = logpfun(x,rho)
f = -abs(x).^rho - log(2) - gammaln(1+1/rho);

function f = ffun(x,rho)
f = rho * sign(x).*abs(x).^(rho-1);
