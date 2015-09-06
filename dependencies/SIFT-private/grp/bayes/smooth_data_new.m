function [B,Theta_diag,C_diag,Alpha_diag,Fit_diag,...
    C_offdiag,Theta_offdiag,Alpha_offdiag,Fit_offdiag phi_t]=smooth_data(C,K,Q,nloop_diag,nloop_offdiag)

%% This function is a wrapper for the smooth_mcmc function, which does the actual smoothing
%% This was written to smoot the diagonal and off-diagonal conectivity time series separately
%% and to save the coefficents and basis functions from each


N=size(C,1);   
M_i=zeros(N,1);
for i=1:N      
    M_i(i)=size(C{i},1);
end


%% Get smoothed time-varying connectivity coefficients

% smooth diagonal connectivities
C_diag=cell(sum(M_i),1);
ind=0;
for i=1:N
    for j=1:M_i(i)
        ind=ind+1;
        C_diag{ind}=squeeze(C{i}(j,j,:));
    end
end    
[ALPHA ALPHA_BAR THETA FIT phi_t]=smooth_mcmc(C_diag,K,Q,nloop_diag);
first=round(nloop_diag/2); last=nloop_diag+1;
Alpha_diag=mean(ALPHA(:,:,first:last),3);
Theta_diag=mean(THETA(:,:,first:last),3);
Fit_diag=phi_t'*Theta_diag*Alpha_diag;

B=cell(N,1);
ind=0;
for i=1:N
    B{i}=cell(M_i(i));
    for j=1:M_i(i)
        ind=ind+1;
        B{i}{j,j}=Alpha_diag(:,ind)';
    end
end


% smooth off-diagonal connectivities
M_offdiag=M_i(1)^2-M_i(1);
for i=2:N
    M_offdiag=M_offdiag+M_i(i)^2-M_i(i);
end
    
C_offdiag=cell(M_offdiag,1);
ind=0;
for i=1:N
    for j1=1:M_i(i)
        for j2=1:M_i(i)
            if(not(j1==j2))
                ind=ind+1;
                C_offdiag{ind}=squeeze(C{i}(j1,j2,:));
            end
        end
    end
end    
[ALPHA ALPHA_BAR THETA FIT phi_t]=smooth_mcmc(C_offdiag,K,Q,nloop_offdiag);
first=round(nloop_offdiag/2); last=nloop_offdiag+1;
Alpha_offdiag=mean(ALPHA(:,:,first:last),3);
Theta_offdiag=mean(THETA(:,:,first:last),3);
Fit_offdiag=phi_t'*Theta_offdiag*Alpha_offdiag;

ind=0;
for i=1:N
    for j1=1:M_i(i)
        for j2=1:M_i(i)
            if not(j1==j2)
                ind=ind+1;
                B{i}{j1,j2}=Alpha_offdiag(:,ind)';
            end
        end
    end
end



