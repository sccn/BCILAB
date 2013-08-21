function [A,W,S,khinds,c,LL,Ltall,gm,alpha,mu,beta,rho] = amica10(x,M,m,maxiter,do_sphere,do_newton,As,cs,Ainit,cinit,kinit)
% function [A,W,S,khinds,c,LL,Ltall,gm,alpha,mu,beta,rho] = amica10(x,[M,m,maxiter,update_rho,mindll,iterwin,do_newton,As,cs])
%
% Performa AMICA algorithm for adaptive mixutre ICA on input data x (of dimension n x N), using M models,
% with adaptive source densities unsing a mixture of m Generalized Gaussians, with fixed or adaptive shape.
%
% Components are shered 
%
% Stopping condition is when avg log likelihood increase is less than mindll (averaged over iterwin iterations)
% and when norm of updates is smaller than minndw.
%
% Inputs:
%
%       x           --  data, channels (n) by time points (N)
%       M           --  number of ICA mixture models (default is 1)
%       m           --  number of source density mixtures (default is 3)
%       maxiter     --  maximum number of iterations to perform (default is 100)
%       do_sphere   --  1 = remove mean and sphere data (default)
%                       2 = only remove mean and normalize channels
%                       0 = do not remove mean or sphere
%       do_newton   --  1 = use Newton update for upmixing matrices (default)
%                       0 = use natural gradient updates
%       As          --  (optional) true generating basis for plotting purposes
%       cs          --  (optional) true model centers for plotting purposes
%
% Outputs:
%
%       A           --  learned mixing matrix with normalized columns, pinv(A) is unmixing matrix ( n x nM )
%       W           --  learned unmixing matrices for sphered data (n x n x M)
%       S           --  Sphering matrix (n x n)
%       khinds      --  column indices of the model components in the returned A matrix (n x M)
%       c           --  learned model centers  ( n x M )
%       LL          --  history of likelihood over iterations ( 1 x min(iter,maxiter) )
%       Ltall       --  log likelihood of time point for each model ( M x N )
%       gm          --  learned ica model mixture proportions ( 1 x M )
%       alpha       --  learned source density mixture proportions ( m x n x M )
%       mu          --  learned source density location parameter ( m x n x M )
%       beta        --  learned source density inverse scale parameter ( m x n x M )
%       rho         --  learned source density shape paramters ( m x n x M )
%
%       A(:,khinds(:,h)) is the mixing matrix associated model h.
%       c(:,h) is the center, or mean, or model h.
%       W(:,:,h)*S is the unmixing matrix associated with model h.
%       Ltall(:,t) is the vector of model likelihoods for time point t
%       gm(h) is the percent of the data explained by model h (associated
%           with the time points where Ltall(h,t) is maximum of Ltall(:,t)
%
%       To plot the likelihood over time points, you can use:
%       
%               >> plot(Ltall');
%
%       or a moving average filtered version of Ltall.
%
%       The probability density of source i in model h (corresponding to
%       component A(:,khinds(i,h)) is:
%
%               p(y) = sum_{j=1}^m alpha(j,i,h) * sqrt(beta(j,i,h)) *
%                   gg(sqrt(beta(j,i,h))*(y-mu(j,i,h)),rho(j,i,h))
%
%       where:
%
%               gg(y,rho(j,i,h)) = (1/(2*GAMMA(1+1/rho(j,i,h)))) *
%                   exp(-abs(y)^rho(j,i,h))
%
%       In EEGLAB, to assign model h to a dataset loaded in the gui, use:
%
%           >> EEG.icaweights = W(:,:,h);
%           >> EEG.icasphere = S;
%           >> eeglab redraw;
%
%       and then possibly save the dataset.
%
% Author: Jason Palmer (jason@sccn.ucsd.edu), Swartz Center for Computational Neuroscience, 2010.
%

comp_thresh = 0.99; % correlation threshold for identifying comps as same
ident_intvl = 100; % iteration interval to check for identical comps

update_gm = 1;
update_alpha = 1;
update_A = 1;
update_mu = 1;
update_beta = 1;
update_rho = 1;
do_reparm = 1;
reparm_mu_c = 0;
fix_init = 1;

if nargin < 5
    do_sphere = 1;
end
if nargin < 4
    maxiter = 200;
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

switch do_sphere
    case 1
        xmn = mean(x,2);
        %S = pinv(sqrtm(cov(x')));
        for i = 1:n
            x(i,:) = x(i,:) - xmn(i);
        end
        [Us,Ss,Vs] = svd(x*x'/N);
        S = Us * diag(1./sqrt(diag(Ss))) * Us';
        x = S*x;
    case 2
        xmn = mean(x,2);
        for i = 1:n
            x(i,:) = x(i,:) - xmn(i);
        end
        S = diag( 1 ./ sqrt(var(x')) );
        x = S*x;
    otherwise
        xmn = mean(x,2);
        for i = 1:n
            x(i,:) = x(i,:) - xmn(i);
        end
        S = eye(n);
end
ldetS = log(abs(det(S)));
if nargin >= 7
    for h = 1:M
        As(:,:,h) = S * As(:,:,h);
        cs(:,h) = S * (cs(:,h) - xmn);
    end
end
lrate0 = 0.1;
lratemax = 1.0;
lnatrate = 0.1;
newt_start_iter = 50;
rho_start_iter = 1;
lratefact = 0.5;

block_size = 3000;


lrate = lrate0;
mindll = 1e-8;
iterwin = 1;

dispinterval = 50;



showLL = 1; % display the log likelihood at each iteration
if nargin >= 7
    plotfig = 1; % we are given the solution bases and centers, As, cs
else
    plotfig = 0;  % can also set this to 1 to plot anyway (without the solution)
end

Mx = max(max(abs(x)));

rho0 = 1.5;
rholrate = 0.1;
rhomin = 1;
rhomax = 2;
if nargin < 6
    do_newton = 1;
end
if nargin >= 7
    plothist = 1;
    nbins = 50;
else
    plothist = 0;
end


% initialize parameters
if nargin >= 9
    A = Ainit;
    num_comps = size(A,2);
    if nargin >= 10
        c = cinit;
    else
        c = zeros(n,M);
    end
    if nargin >= 11
        khinds = kinit;
    else
        khinds = reshape(1:(n*M),n,M);
    end
else
    num_comps = n*M;
    for h = 1:M
        A(:,(h-1)*n+1:h*n) = eye(n) + 0.05*(0.5-rand(n,n));
        c(:,h) = randn(n,1);
        khinds(:,h) = ((h-1)*n+1):(h*n);
    end
    for k = 1:num_comps
        A(:,k) = A(:,k) / norm(A(:,k));
    end
end

gm = (1/M) * ones(M,1);
alpha = (1/m) * ones(m,n*M);

if fix_init
    mu = ( (1:m)'-1-(m-1)/2 ) * ones(1,n*M);
    beta = ones(m,n*M);
else
    if m > 1
        mu = 0.1 * randn(m,n*M);
    else
        mu = zeros(m,n*M);
    end
    beta = ones(m,n*M) + 0.1 * rand(m,n*M);
end
rho = rho0 * ones(m,n*M);

% allocate variables
b = zeros(n,block_size,M);
g = zeros(n,block_size);

y = zeros(n,block_size,m,M);
fp = zeros(block_size);
Q = zeros(m,block_size);

if nargin >= 7
    Ltall = zeros(M,N);
end
Lt = zeros(M,block_size);
v = zeros(M,block_size);
z = zeros(n,block_size,m,M);

num_blocks = max(1,floor(N / block_size));

% iteratively optimize parameters using EM
for iter = 1:maxiter
    for h = 1:M
        W(:,:,h) = pinv(A(:,khinds(:,h)));
        ldet(h) = log(abs(det(W(:,:,h))));
    end
        
    % Loop over data blocks and accumulate updates
    LL(iter) = 0;
    dalpha_numer = zeros(m,num_comps);
    dalpha_denom = zeros(m,num_comps);
    dmu_numer = zeros(m,num_comps);
    dmu_denom = zeros(m,num_comps);
    dbeta_numer = zeros(m,num_comps);
    dbeta_denom = zeros(m,num_comps);
    drho_numer = zeros(m,num_comps);
    drho_denom = zeros(m,num_comps);

    vsumsum = zeros(1,M);

    dbaralpha_numer = zeros(m,n,M);
    dbaralpha_denom = zeros(m,n,M);
    dlambda_numer = zeros(m,n,M);
    dlambda_denom = zeros(m,n,M);
    dsigma2_numer = zeros(n,M);
    dsigma2_denom = zeros(n,M);
    dkappa_numer = zeros(m,n,M);
    dkappa_denom = zeros(m,n,M);
    Phi = zeros(n,n,M);
    cnew = zeros(n,M);
    
    blksize = block_size;
    for blk = 1:num_blocks
        xstrt = (blk-1)*block_size + 1;
        if blk == num_blocks
            blksize = size(x,2) - xstrt + 1;
        end
        xstp = xstrt + blksize - 1;
                
        % get y, z, Q, Lt
        for h = 1:M
            Lt(h,1:blksize) = log(gm(h)) + ldet(h) + ldetS;
            for i = 1:n
                b(i,1:blksize,h) = W(i,:,h) * x(:,xstrt:xstp) - W(i,:,h)*c(:,h);
                for j = 1:m
                    y(i,1:blksize,j,h) = sqrt(beta(j,khinds(i,h))) * ( b(i,1:blksize,h) - mu(j,khinds(i,h)) );
                    Q(j,1:blksize) = log(alpha(j,khinds(i,h))) + 0.5*log(beta(j,khinds(i,h))) + ...
                        feval('logpfun',y(i,1:blksize,j,h),rho(j,khinds(i,h)));
                end
                Qmax(1:blksize) = max(Q(:,1:blksize),[],1);
                qtmp(1:blksize) = zeros(1,blksize);
                for j = 1:m
                    qtmp(1:blksize) = qtmp(1:blksize) + exp(Q(j,1:blksize) - Qmax(1:blksize));
                end
                tmpvec(1:blksize) = Qmax(1:blksize) + log(qtmp(1:blksize));
                Lt(h,1:blksize) = Lt(h,1:blksize) + tmpvec(1:blksize);                                    
                for j = 1:m
                    z(i,1:blksize,j,h) = 1 ./ exp(tmpvec(1:blksize)-Q(j,1:blksize));
                end
            end
        end
        
        Ltall(:,xstrt:xstp) = Lt(:,1:blksize);
        
        Ltmax(1:blksize) = max(Lt(:,1:blksize),[],1);
        
        vtmp(1:blksize) = zeros(1,blksize);
        for h = 1:M
            vtmp(1:blksize) = vtmp(1:blksize) + exp(Lt(h,1:blksize) - Ltmax(1:blksize));
        end
        P(1:blksize) = Ltmax(1:blksize) + log(vtmp(1:blksize));
        LL(iter) = LL(iter) + sum(P(1:blksize));
        for h = 1:M
            v(h,1:blksize) = 1 ./ exp(P(1:blksize) - Lt(h,1:blksize));
        end
        
        % get v, u, gm, and alpha
        for h = 1:M
            vsum = sum(v(h,1:blksize));
            vsumsum(h) = vsumsum(h) + vsum;
            
            cnew(:,h) = cnew(:,h) + x(:,xstrt:xstp) * v(h,1:blksize)';
            g(:,1:blksize) = zeros(n,blksize);
            
            for i = 1:n
                dsigma2_numer(i,h) = dsigma2_numer(i,h) + sum(v(h,1:blksize).*b(i,1:blksize,h).^2);
                dsigma2_denom(i,h) = dsigma2_denom(i,h) + vsum;
                
                for j = 1:m
                    u(1:blksize) = v(h,1:blksize) .* z(i,1:blksize,j,h);
                    usum = sum(u(1:blksize));
                                        
                    fp(1:blksize) = feval('ffun',y(i,1:blksize,j,h),rho(j,khinds(i,h)));
                    ufp(1:blksize) = u(1:blksize) .* fp(1:blksize);
                    
                    g(i,1:blksize) = g(i,1:blksize) + sqrt(beta(j,khinds(i,h))) * ufp(1:blksize);
                                                            
                    dkappa_numer(j,i,h) = dkappa_numer(j,i,h) + beta(j,khinds(i,h)) * sum(ufp(1:blksize).*fp(1:blksize));
                    dkappa_denom(j,i,h) = dkappa_denom(j,i,h) + usum;
                    
                    dlambda_numer(j,i,h) = dlambda_numer(j,i,h) + sum(u(1:blksize).*(fp(1:blksize).*y(i,1:blksize,j,h)-1).^2);
                    dlambda_denom(j,i,h) = dlambda_denom(j,i,h) + usum;

                    dalpha_numer(j,khinds(i,h)) = dalpha_numer(j,khinds(i,h)) + usum;
                    dalpha_denom(j,khinds(i,h)) = dalpha_denom(j,khinds(i,h)) + vsum;

                    dbaralpha_numer(j,i,h) = dbaralpha_numer(j,i,h) + usum;
                    dbaralpha_denom(j,i,h) = dbaralpha_denom(j,i,h) + vsum;
                    
                    if rho(j,khinds(i,h)) <= 2
                        dmu_numer(j,khinds(i,h)) = dmu_numer(j,khinds(i,h)) + sum(ufp(1:blksize));
                        dmu_denom(j,khinds(i,h)) = dmu_denom(j,khinds(i,h)) + ...
                            sqrt(beta(j,khinds(i,h))) * sum(ufp(1:blksize)./y(i,1:blksize,j,h));

                        dbeta_numer(j,khinds(i,h)) = dbeta_numer(j,khinds(i,h)) + usum;
                        dbeta_denom(j,khinds(i,h)) = dbeta_denom(j,khinds(i,h)) + sum(ufp(1:blksize).*y(i,1:blksize,j,h));
                    else
                        dmu_numer(j,khinds(i,h)) = dmu_numer(j,khinds(i,h)) + sum(ufp(1:blksize));
                        dmu_denom(j,khinds(i,h)) = dmu_denom(j,khinds(i,h)) + ...
                            sqrt(beta(j,khinds(i,h))) * sum(ufp(1:blksize).*fp(1:blksize));

                        %db = ( rho(j,inds(i,h)) * sum(z(i,:,j,h).*abs(y(i,:,j,h)).^rho(j,inds(i,h))) )^(-2 / rho(j,inds(i,h)));
                        dbeta_denom(j,khinds(i,h)) = dbeta_denom(j,khinds(i,h)) + ...
                            rho(j,khinds(i,h)) * sum(u(1:blksize).*abs(y(i,1:blksize,j,h)).^rho(j,khinds(i,h)));
                    end
                    
                    if update_rho
                        ytmp0(1:blksize) = abs(y(i,1:blksize,j,h)).^rho(j,khinds(i,h));
                        ytmp(1:blksize) = ytmp0(1:blksize) .* log(ytmp0(1:blksize));
                        %ytmp(ytmp0(1:blksize)<=1e-8) = 0;
                        %drho_numer(j,khinds(i,h)) = drho_numer(j,khinds(i,h)) + sum(u(1:blksize).*ytmp(1:blksize).*log(ytmp(1:blksize)));
                        drho_numer(j,khinds(i,h)) = drho_numer(j,khinds(i,h)) + sum(u(1:blksize) .* ytmp(1:blksize));
                        drho_denom(j,khinds(i,h)) = drho_denom(j,khinds(i,h)) + usum;

                        %dr2 = 1 - rho(j,inds(i,h)) * dr / psi(1+1/rho(j,inds(i,h)));
                    end
                end
            end
            Phi(:,:,h) = Phi(:,:,h) + g(:,1:blksize) * b(:,1:blksize,h)';
            
        end
        
    end

    LL(iter) = LL(iter)/(n*N);

    if showLL
        disp(['iter '  int2str(iter) ' lrate = ' num2str(lrate) ' LL = ' num2str(LL(iter))]);
    end
    if iter > 1 && LL(iter) < LL(iter-1)
        disp('Likelihood decreasing!');
        lrate = lrate * lratefact;
        %lratemax = lratemax * lratefact;
        %rholrate = rholrate * lratefact;
    end

    % Do parameter updates
    if update_gm
        gm = vsumsum / N;
    end
    if update_alpha
        alpha = dalpha_numer ./ dalpha_denom;
    end
    baralpha = dbaralpha_numer ./ dbaralpha_denom;    
    sigma2 = dsigma2_numer ./ dsigma2_denom;    
    kappa = zeros(n,h);
    lambda = zeros(n,h);
    for h = 1:M
        cnew(:,h) = cnew(:,h) / vsumsum(h);
        for i = 1:n
            for j = 1:m
                dk = dkappa_numer(j,i,h) ./ dkappa_denom(j,i,h);
                kappa(i,h) = kappa(i,h) + baralpha(j,i,h) * dk;
                lambda(i,h) = lambda(i,h) + baralpha(j,i,h) * ...
                    ( dlambda_numer(j,i,h) ./ dlambda_denom(j,i,h) + dk * mu(j,khinds(i,h))^2 );
            end
        end
        Phi(:,:,h) = Phi(:,:,h) / vsumsum(h);
    end
    
    if update_mu
        mu = mu + dmu_numer ./ dmu_denom;
    end
    if update_beta
        beta = beta .* dbeta_numer ./ dbeta_denom;
    end
    if update_rho && iter >= rho_start_iter
        rho = rho + rholrate * ( 1 - (rho./psi(1+1./rho)) .* drho_numer ./ drho_denom );
        %rho = rho + rholrate * ( psi(1+1./rho)./rho - drho_numer ./ drho_denom );
        rho = min(rho,rhomax);
        rho = max(rho,rhomin);
    end
    
    % update A
    no_newt = 0;
    for h = 1:M
        dAtmp = eye(n) - Phi(:,:,h);
        
        bflag = 0;
        for i = 1:n
            for k = 1:n
                if i == k
                    B(i,i) = dAtmp(i,i) / lambda(i,h);
                else
                    denom = kappa(i,h)*kappa(k,h)*sigma2(i,h)*sigma2(k,h) - 1;
                    if denom > 0
                        B(i,k) = (-kappa(k,h) * sigma2(i,h) * dAtmp(i,k) + dAtmp(k,i)) / denom;
                    else
                        bflag = 1;
                        no_newt = 1;
                    end
                end
            end
        end
        if bflag == 0 && do_newton == 1 && iter >= newt_start_iter
            dA(:,:,h) = A(:,khinds(:,h)) * B;
        else
            if bflag == 1 && iter >= newt_start_iter
                disp('Hessian non positive definite! Using natural gradient');
            end
            dA(:,:,h) = -A(:,khinds(:,h)) * dAtmp;
        end        
    end
        
    dAk = zeros(n,num_comps);
    zeta = zeros(num_comps);
    for h = 1:M
        for i = 1:n
            dAk(:,khinds(i,h)) = dAk(:,khinds(i,h)) + gm(h)*dA(:,i,h);
            zeta(khinds(i,h)) = zeta(khinds(i,h)) + gm(h);
        end
    end
    for k = 1:num_comps
        dAk(:,k) = dAk(:,k) / zeta(k);
    end
        
    cmps = unique(khinds(:));
    if update_A
        if iter >= newt_start_iter && do_newton && (no_newt==0)
            lrate = min(lratemax,lrate + min(0.1,lrate));
        else
            lrate = min(lnatrate,lrate + min(0.1,lrate));
        end
        A = A + lrate * dAk;
    
        % reparameterize
        % normalize components
        if do_reparm
            for k = 1:length(cmps)
                tau = norm(A(:,cmps(k)));
                A(:,cmps(k)) = A(:,cmps(k)) / tau;
                mu(:,cmps(k)) = mu(:,cmps(k)) * tau;
                beta(:,cmps(k)) = beta(:,cmps(k)) / tau^2;
                if reparm_mu_c
                    %dmu(cmps(k)) = sum(alpha(:,cmps(k)).*mu(:,cmps(k)));
                    % get a model and source index associated with this component
                    tmp = find(khinds(:)==cmps(k),1);
                    hk = ceil(tmp / n);
                    ik = mod(tmp,n);
                    if ik == 0
                        ik = n;
                    end
                    dmu = W(ik,:,hk) * (cnew(:,hk) - c(:,hk));
                    for j = 1:m
                        mu(j,cmps(k)) = mu(j,cmps(k)) - dmu;
                    end
                end
            end
            c = cnew;
        end
    end
    
    % identify components
    if mod(iter,ident_intvl) == 0
        for h = 1:M
            for hh = (h+1):M
                for i = 1:n
                    for ii = 1:n
                        if abs(A(:,khinds(i,h))'*A(:,khinds(ii,hh))) >= comp_thresh && khinds(i,h) ~= khinds(ii,hh)
                            if ~any(khinds(:,hh)==khinds(i,h)) % check if the similar component is already identified with another comp in model
                                disp(['identifying comp ' int2str(khinds(ii,hh)) ' with comp ' int2str(khinds(i,h))]);
                                khinds(khinds(:)==khinds(ii,hh)) = khinds(i,h);
                            end
                        end
                    end
                end
            end
        end
    end
                            
        
    % Plot
    if plotfig && mod(iter,dispinterval) == 0
        figure(1), clf, axis(Mx*[-1 1 -1 1]);hold on;
        plot(x(1,:),x(2,:),'*');

        for h = 1:M
            for i = 1:n
                if nargin >= 7
                    plot([cs(1,h) cs(1,h)+As(1,i,h)],[cs(2,h) cs(2,h)+As(2,i,h)],'r');
                end
                plot([c(1,h) c(1,h)+A(1,khinds(i,h))],[c(2,h) c(2,h)+A(2,khinds(i,h))],'g');
            end
        end
        pause(0.5);
    end
    
    if plothist && mod(iter,dispinterval) == 0
        [mval,mind] = max(Ltall,[],1);
        for h = 1:M
            indsh = find(mind == h);
            
            figure(2); 
            for i = 1:n
                ball = W(i,:,h) * x - W(i,:,h)*c(:,h);
                [hst,bn] = hist(ball(indsh),nbins);

                dbn = abs(bn(2)-bn(1));
                hst = hst / (sum(hst)*dbn);

                ph = zeros(1,length(bn));

                subplot(M,n,n*(h-1)+i); hold off; plot(0);

                bmax = max(abs(bn(1)),abs(bn(end)));
                axis([-bmax-2*dbn bmax+2*dbn 0 max(hst)+0.1]); hold on;
                for j = 1:m
                    dh = alpha(j,khinds(i,h)) * sqrt(beta(j,khinds(i,h))) * exp(-abs(sqrt(beta(j,khinds(i,h)))*(bn-mu(j,khinds(i,h)))).^rho(j,khinds(i,h))) / (2*gamma(1+1/rho(j,khinds(i,h))));
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

A = pinv(S)*A;
c = pinv(S)*c;
for k = 1:num_comps
    tau = norm(A(:,k));
    A(:,k) = A(:,k) / tau;
    mu(:,k) = mu(:,k) * tau;
    beta(:,k) = beta(:,k) / tau^2;
end




function f = pfun(x,rho)
f = (1 / (2*GAMMA(1+1/rho))) * exp( - abs(x).^rho );

function f = logpfun(x,rho)
f = -abs(x).^rho - log(2) - gammaln(1+1/rho);

function f = ffun(x,rho)
f = rho * sign(x).*abs(x).^(rho-1);
