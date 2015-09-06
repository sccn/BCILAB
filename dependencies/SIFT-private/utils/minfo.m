function [MI,vMI,MInrm,Hu,vHu,Hx,vHx] = minfo(u,x,nbins,nlags,verb)
% function MI = minfo(u,x,nbins,nlags) : Calculate the pairwise mutual
% information between rows of u, or between rows of u and rows of x if
% second argument is a matrix. If nlags > 0, calculate for lags
% -nlags:nlags. The output is a matrix of size [nu,nx,2*nlags+1] of
% pairwise mutual information values, with zero lag at MI(:,:,nlags+1).
% Uses nbins bins for pdf of u(i,:), equally spaced between min(u(i,:))
% and max(u(i,:)), similarly for pdfs of x rows.
%
% Inputs:
%           u       Matrix (nu by Nu) of nu time series.
%           x       Optional second input matrix (nx by Nx) of nx time
%                   series. Nx must equal Nu, i.e. all time series must
%                   have the same length. Enter 0 or [] to bypass in order
%                   to set nbins or nlags.
%           nbins   Number of bins to use in computing pdfs. Default is
%                   round(3*log2(1+Nu/10)). Use [] or 0 to get default.
%           nlags   Number of lags (in addition to zero lag) to compute
%                   mutual information (before and after zero) lag.
%                   Default is 0.
%           verb    Flag for verbose output. Def == 0.
%
% Outputs:
%           MI      If no second input x, or x==0, then MI is an
%                   (nu by nu by 2*nlags+1) matrix of pairwise mutual
%                   information between rows of u at lags -nlags to nlags.
%                   If both inputs u and x are present, then MI is an
%                   (nu by nx by 2*nlags+1) matrix of pairwise mutual
%                   information between rows of u and rows of x. For
%                   example, MI(i,j,nlags+1+g) is the mutual information
%                   between u(i,1+g:N+g) and u(j,1:N), or between
%                   u(i,1+g:N+g) and x(j,1:N). So a peak at g+nlags+1 in
%                   MI(i,j,:) means i lags j by g time points (i leads j by
%                   -g time points.)
%           MInrm   Matrix (nu by nu by 2*nlags+1) of "normalized" pairwise
%                   mutual information of rows of u, relative to the
%                   average of the marginal (non-differential) bin
%                   entropies. If x argument is present, normalization is
%                   always with respect to u(i,:).
%           Hu      Vector nu by 1 of differential entropies of rows of u.
%           Hx      Vector nu by 1 of differential entropies of rows of x.
%           vXX     Variance in estimate of XX
%           
% Author: Jason Palmer, SCCN/INC, 2009

varrng = 0;
verb = 0;
secndord = 0;
[nu,Nu] = size(u);
if nargin < 3 || isempty(nbins) || nbins == 0
    nbins = round(3*log2(1+Nu/10));
end
if nargin < 4
    nlags = 0;
end
if nargin >= 2
    [nx,Nx] = size(x);
    if Nx == 1 || Nx == 0
        MI = zeros(nu,nu,2*nlags+1); 
        MInrm = zeros(nu,nu,2*nlags+1);
        dox = 0;
    else        
        if Nx ~= Nu
            error('input time series must be same length');
        end
        MI = zeros(nu,nx,2*nlags+1);
        MInrm = zeros(nu,nx,2*nlags+1);        
        dox = 1;
    end
else
    MI = zeros(nu,nu,2*nlags+1);
    MInrm = zeros(nu,nu,2*nlags+1);
    dox = 0;
end

% bin the time series
if verb
    disp('binning the time series ...'); pause(0.1);
end
Hu = zeros(nu,2*nlags+1);
deltau = zeros(nu,1);
for i = 1:nu
    if varrng
        um = mean(u(i,:));
        us = std(u(i,:));
        umax = min(max(u(i,:)),um + 5*us);
        umin = max(min(u(i,:)),um - 5*us);
    else
        umax = max(u(i,:));
        umin = min(u(i,:));
    end
    deltau(i) = (umax-umin)/nbins;
    u(i,:) = 1 + round((nbins - 1) * (u(i,:) - umin) / (umax - umin));
    u(i,:) = min(nbins,u(i,:));
    u(i,:) = max(1,u(i,:));
    
    if nlags == 0
        pmfr = diff([0 find(diff(sort(u(i,:)))) Nu])/Nu;
        Hu(i) = -sum(pmfr.*log(pmfr));
        vHu(i) = (sum(pmfr.*(log(pmfr).^2)) - Hu(i)^2) / Nu;
        if secndord
            Hu(i) = Hu(i) + (nbins-1)/(2*Nu) - (1/12/Nu^2)*(1-sum(1./pmfr));
        else
            Hu(i) = Hu(i) + (nbins-1)/(2*Nu);
        end
    else
        for g = 0:nlags
            pmfr = diff([0 find(diff(sort(u(i,1+g:Nu)))) (Nu-g)])/(Nu-g);    
            Hu(i,nlags+1+g) = -sum(pmfr.*log(pmfr));
            vHu(i,nlags+1+g) = (sum(pmfr.*(log(pmfr).^2)) - Hu(i,nlags+1+g)^2) / Nu;
            if secndord
                Hu(i,nlags+1+g) = Hu(i,nlags+1+g) + (nbins-1)/(2*Nu) - (1/12/Nu^2)*(1-sum(1./pmfr));
            else
                Hu(i,nlags+1+g) = Hu(i,nlags+1+g) + (nbins-1)/(2*Nu);
            end
            if g > 0
                pmfr = diff([0 find(diff(sort(u(i,1:Nu-g)))) (Nu-g)])/(Nu-g);    
                Hu(i,nlags+1-g) = -sum(pmfr.*log(pmfr));
                vHu(i,nlags+1-g) = (sum(pmfr.*(log(pmfr).^2)) - Hu(i,nlags+1-g)^2) / Nu;
                if secndord
                    Hu(i,nlags+1-g) = Hu(i,nlags+1-g) + (nbins-1)/(2*Nu) - (1/12/Nu^2)*(1-sum(1./pmfr));
                else
                    Hu(i,nlags+1-g) = Hu(i,nlags+1-g) + (nbins-1)/(2*Nu);
                end
            end
        end
    end
end
if dox == 1
    Hx = zeros(nx,2*nlags+1);
    deltax = zeros(nx,1);
    for i = 1:nx
        xmax = max(x(i,:));
        xmin = min(x(i,:));
        deltax(i) = (xmax-xmin)/nbins;
        x(i,:) = 1 + round((nbins - 1) * (x(i,:) - xmin) / (xmax - xmin));

        if nlags == 0
            pmfr = diff([0 find(diff(sort(x(i,:)))) Nx])/Nx;                        
            Hx(i) = -sum(pmfr.*log(pmfr));
            vHx(i) = (sum(pmfr.*(log(pmfr).^2)) - Hx(i)^2) / Nx;
            if secndord
                Hx(i) = Hx(i) + (nbins-1)/(2*Nx) - (1/12/Nx^2)*(1-sum(1./pmfr));
            else
                Hx(i) = Hx(i) + (nbins-1)/(2*Nx);
            end
        else
            for g = 0:nlags
                pmfr = diff([0 find(diff(sort(x(i,1+g:Nx)))) (Nx-g)])/(Nx-g);
                Hx(i,nlags+1+g) = -sum(pmfr.*log(pmfr));
                vHx(i,nlags+1+g) = (sum(pmfr.*(log(pmfr).^2)) - Hx(i,nlags+1+g)^2) / Nx;
                if secndord
                    Hx(i,nlags+1+g) = Hx(i,nlags+1+g) + (nbins-1)/(2*Nx) - (1/12/Nx^2)*(1-sum(1./pmfr));
                else
                    Hx(i,nlags+1+g) = Hx(i,nlags+1+g) + (nbins-1)/(2*Nx);
                end
                if g > 0
                    pmfr = diff([0 find(diff(sort(x(i,1:Nx-g)))) (Nx-g)])/(Nx-g);    
                    Hx(i,nlags+1-g) = -sum(pmfr.*log(pmfr));
                    vHx(i,nlags+1-g) = (sum(pmfr.*(log(pmfr).^2)) - Hx(i,nlags+1-g)^2) / Nx;
                    if secndord
                        Hx(i,nlags+1-g) = Hx(i,nlags+1-g) + (nbins-1)/(2*Nx) - (1/12/Nx^2)*(1-sum(1./pmfr));
                    else
                        Hx(i,nlags+1-g) = Hx(i,nlags+1-g) + (nbins-1)/(2*Nx);
                    end
                end
            end
        end
    end
else
    Hx = 0;
    vHx = 0;
end

% get pairwise histograms at each lag
if verb
    disp('getting pairwise mutual information ...'); pause(0.1);
end
if dox == 0
    for i = 1:nu
        if i == 1
            if verb
                disp(['doing ' int2str(i) ' ...']); pause(0.01);
            end
            tic;
        else
            t1 = toc; T = t1 * (nu-i+1);
            if verb
                disp(['doing ' int2str(i) ' ... time rem: ' num2str(T/60) ' m']); pause(0.01);
            end
            tic;
        end
        for j = 1:nu
            if nlags == 0 % faster if we don't need to compute lag index                
                if i == j
                    MI(i,i) = Hu(i);
                    MInrm(i,i) = 1;
                else
                    pmf2r = diff([0 find(diff(sort(u(i,:)+nbins*(u(j,:)-1)))) Nu])/Nu;            
                    H2 = -sum(pmf2r.*log(pmf2r));
                    MI(i,j) = Hu(i) + Hu(j) - H2 - (nbins^2-1)/2/Nu;
                    MI(j,i) = MI(i,j);

                    vMI(i,j) = vHu(i) + vHu(j) + (sum(pmf2r.*(log(pmf2r).^2)) - H2^2) / Nu;
                    vMI(j,i) = vMI(i,j);

                    MInrm(i,j) = MI(i,j) / Hu(i);
                    MInrm(j,i) = MI(j,i) / Hu(j);
                end                
            else
                % do zero lag
                if i == j
                    MI(i,i,nlags+1) = Hu(i,nlags+1);
                    MInrm(i,i,nlags+1) = 1;
                else                    
                    pmf2r = diff([0 find(diff(sort(u(i,1:Nu)+nbins*(u(j,1:Nu)-1)))) Nu])/Nu;
                    H2 = -sum(pmf2r.*log(pmf2r));
                    MI(i,j,nlags+1) = Hu(i,nlags+1) + Hu(j,nlags+1) - H2 - (nbins^2-1)/2/Nu;
                    MI(j,i,nlags+1) = MI(i,j,nlags+1);

                    vMI(i,j,nlags+1) = vHu(i,nlags+1) + vHu(j,nlags+1) + (sum(pmf2r.*(log(pmf2r).^2)) - H2^2) / Nu;
                    vMI(j,i,nlags+1) = vMI(i,j,nlags+1);

                    MInrm(i,j,nlags+1) = MI(i,j,nlags+1) / Hu(i,nlags+1);
                    MInrm(j,i,nlags+1) = MI(j,i,nlags+1) / Hu(j,nlags+1);
                end

                for g = 1:nlags
                    pmf2r = diff([0 find(diff(sort(u(i,1+g:Nu)+nbins*(u(j,1:(Nu-g))-1)))) Nu-g])/(Nu-g);
                    H2 = -sum(pmf2r.*log(pmf2r));
                    MI(i,j,nlags+1+g) = Hu(i,nlags+1+g) + Hu(j,nlags+1-g) - H2 - (nbins^2-1)/2/Nu;

                    vMI(i,j,nlags+1+g) = vHu(i,nlags+1+g) + vHu(j,nlags+1-g) + (sum(pmf2r.*(log(pmf2r).^2)) - H2^2) / Nu;
                      
                    MInrm(i,j,nlags+1+g) = MI(i,j,nlags+1+g) / Hu(i,nlags+1+g);
                    if i == j
                        pmf2r = diff([0 find(diff(sort(u(i,1:Nu-g)+nbins*(u(j,1+g:Nu)-1)))) Nu-g])/(Nu-g);
                        H2 = -sum(pmf2r.*log(pmf2r));
                        MI(i,j,nlags+1-g) = Hu(i,nlags+1+g) + Hu(j,nlags+1-g) - H2 - (nbins^2-1)/2/Nu;

                        vMI(i,j,nlags+1-g) = vHu(i,nlags+1+g) + vHu(j,nlags+1-g) + (sum(pmf2r.*(log(pmf2r).^2)) - H2^2) / Nu;
                        
                        MInrm(i,j,nlags+1-g) = MI(i,j,nlags+1-g) / Hu(i,nlags+1-g);
                    else
                        MI(j,i,nlags+1-g) = MI(i,j,nlags+1+g);
                        vMI(j,i,nlags+1-g) = vMI(i,j,nlags+1+g);
                        
                        MInrm(j,i,nlags+1-g) = MI(j,i,nlags+1-g) / Hu(j,nlags+1-g);
                    end
                end                
            end        
        end
    end    
else
    for i = 1:nu
        if i == 1
            if verb
                disp(['doing ' int2str(i) ' ...']); pause(0.01);
            end
            tic;
        else
            t1 = toc; t0 = t1/(nu-i+2); T = t0 * (nu-i+2)*(nu-i+1)/2;
            if verb
                disp(['doing ' int2str(i) ' ... time rem: ' num2str(T/60) ' m']); pause(0.01);
            end
        end
        for j = 1:nx
            if nlags == 0 % faster if we don't need to compute lag index
                pmf2r = diff([0 find(diff(sort(u(i,:)+nbins*(x(j,:)-1)))) Nu])/Nu;            
                H2 = -sum(pmf2r.*log(pmf2r));
                MI(i,j) = Hu(i) + Hx(j) - H2 - (nbins^2-1)/2/Nu;

                vMI(i,j) = vHu(i) + vHx(j) + (sum(pmf2r.*(log(pmf2r).^2)) - H2^2) / Nu;

                MInrm(i,j) = MI(i,j) / Hu(i);
            else            
                for g = 0:nlags
                    pmf2r = diff([0 find(diff(sort(u(i,1+g:Nu)+nbins*(x(j,1:(Nu-g))-1)))) Nu-g])/(Nu-g);     
                    H2 = -sum(pmf2r.*log(pmf2r));
                    MI(i,j,nlags+1+g) = Hu(i,nlags+1+g) + Hx(j,nlags+1-g) - H2 - (nbins^2-1)/2/Nu;
                    
                    vMI(i,j,nlags+1+g) = vHu(i,nlags+1+g) + vHx(j,nlags+1-g) + (sum(pmf2r.*(log(pmf2r).^2)) - H2^2) / Nu;
                    
                    MInrm(i,j,nlags+1+g) = MI(i,j,nlags+1+g) / Hu(i,nlags+1+g);
                
                    if g > 0                
                        pmf2r = diff([0 find(diff(sort(u(i,1:(Nu-g))+nbins*(x(j,1+g:Nu)-1)))) Nu-g])/(Nu-g);        
                        H2 = -sum(pmf2r.*log(pmf2r));
                        MI(i,j,nlags+1-g) = Hu(i,nlags+1-g) + Hx(j,nlags+1+g) - H2 - (nbins^2-1)/2/Nu;

                        vMI(i,j,nlags+1-g) = vHu(i,nlags+1-g) + vHx(j,nlags+1+g) + (sum(pmf2r.*(log(pmf2r).^2)) - H2^2) / Nu;
                        
                        MInrm(i,j,nlags+1-g) = MI(i,j,nlags+1-g) / Hu(i,nlags+1-g);
                    end
                end
            end        
        end        
    end
end
Hu = Hu(:,nlags+1);
for i = 1:nu
    Hu(i) = Hu(i) + log(deltau(i));    
end
if dox == 1
    Hx = Hx(:,nlags+1);
    for i = 1:nx
        Hx(i) = Hx(i) + log(deltax(i));
    end
end
