%% Read in data
clear all; clc;
load cluster_centroids_and_mappings

%subj_keep=[1:12 14:15 17:22 24];    %subjects to keep
subj_keep=[1:24];
N=length(subj_keep);                 %number of subjects
S=cell(N,1);                    %cell array for spatial location data
C=cell(N,1);                    %cell array for time-varying connectivity data
T=80; time=1:T;                 %number of time points in eeg time series
                             
dataPath = '/Users/timmullen/Documents/WORK/Wes/MutiSubjectBayesian/Data/';
%% The following block of code is my attempt to simplify the data by
% removing on of the symmetric dipoles for each pair and removing sources
% which were put in the same cluster.   
ind_keep=cell(N,1);
ind_sym=cell(N,1);
for i=1:N                       %find & remove symmetric components & repeated clusterings
    ii=subj_keep(i);
    x=strcat(dataPath,subjectname{ii},'_theta_Wrong_Correct.mat');
    load(x)
    clust_i=ones(1,length(clustcps{1}{ii}));
    ind_i=clustcps{1}{ii}; % Component numbers for subject i
    for k=2:12
        clust_i=[clust_i k*ones(1,length(clustcps{k}{ii}))];
        ind_i=[ind_i clustcps{k}{ii}];
    end    
    tmp=size(Conn,1);
    sources_i=sources.model(1).component;
    for j=2:tmp
        sources_i=[sources_i sources.model(j).component];
    end
    ind_sym{i}=zeros(tmp,1);
    ind_keep{i}=zeros(tmp,1);
    for j=1:tmp
        if sum(sources_i(j)==ind_i)>0
            ind_keep{i}(j)=1;
        end
        if(size(sources(1).model(j).posxyz,1)>1)
            if not(sources(1).model(j).posxyz(2,1)==0)
                ind_sym{i}(j)=1;
            end
        end
    end
    ind_rep=find(ind_keep{i}==1);
    for j=2:sum(ind_keep{i})
        if clust_i(j)==clust_i(j-1)
            ind_keep{i}(find(sources_i==ind_i(j)))=0;
        end
    end
end


%% read in data and save in cell array S of spatial coordinates for input into MCMC algorithm
M_i=zeros(N,1);
for i=1:N                      
    ii=subj_keep(i);
    x=strcat(dataPath,subjectname{ii},'_theta_Wrong_Correct.mat');
    load(x)
    C{i}=log(Conn(find(ind_keep{i}==1),find(ind_keep{i}==1),:)+1);
    M_i(i)=size(C{i},1);
    S{i}=zeros(M_i(i),3);
    for j=1:M_i(i)
        % only save data for which ind_keep==1 and handle symmetric dipoles
        ind_i=find(ind_keep{i}==1); 
        if ind_sym{i}(ind_i(j))==0
            S{i}(j,:)=sources(1).model(ind_i(j)).posxyz(1,:);
        end
        if and(ind_sym{i}(ind_i(j))==1,sources(1).model(ind_i(j)).posxyz(1,2)<0)
            S{i}(j,:)=sources(1).model(ind_i(j)).posxyz(1,:);
        end        
        if and(ind_sym{i}(ind_i(j))==1,sources(1).model(ind_i(j)).posxyz(1,2)>0)
            S{i}(j,:)=sources(1).model(ind_i(j)).posxyz(2,:);
        end        
    end
end


%% Smooth time-varying connectivity and obtain cell array B of connectivity coefficients
K=8;    %number of basis functions to use for fitting connectivity data
Q=5;    %number of FPCs to keep in fitting 

% get smoothing cofficients for diagonal and off-diagonal connectivities
[fit_distrib, B,Theta_diag,C_diag,Alpha_diag,Fit_diag,...
    C_offdiag,Theta_offdiag,Alpha_offdiag,Fit_offdiag phi_t]=smooth_data_new(C,K,Q,10,10); % return basis coeffients corresponding to time-varying connectivity

% rotate and scale the basis functions for diagonal connectivities so that the resulting
% smoothing coefficents are orthonormal. This is because the mean 
% posterior estimates from smooth_mcmc are not necessarily orthonormal
Alpha_diag_mean=mean(Alpha_diag,2);  % mean of posterior estimates of smoothing coefficients
S_Alpha_diag=(Alpha_diag-repmat(Alpha_diag_mean,1,size(Alpha_diag,2)))*...
    (Alpha_diag-repmat(Alpha_diag_mean,1,size(Alpha_diag,2)))'; % sample covariance of coefs
Sigma_psi_diag=Theta_diag*S_Alpha_diag*Theta_diag';  
Sigma_psi_diag=.5*(Sigma_psi_diag+Sigma_psi_diag');   
[U D]=eig(Sigma_psi_diag);
Theta_diag_new=U(:,(K-Q+1):K)*D((K-Q+1):K,(K-Q+1):K)^.5;
Alpha_diag_new=regress(squeeze(C{1}(1,1,:)),phi_t'*Theta_diag_new);
for i=1:N
    for j=1:M_i(i)
        Alpha_diag_new=[Alpha_diag_new...
            regress(squeeze(C{i}(j,j,:)),phi_t'*Theta_diag_new)];
    end
end
Alpha_diag_new=Alpha_diag_new(:,2:size(Alpha_diag_new,2));
Theta_diag=Theta_diag_new*cov(Alpha_diag_new')^.5;
Alpha_diag=cov(Alpha_diag_new')^-.5*Alpha_diag_new;
Fit_diag=phi_t'*Theta_diag*Alpha_diag;

% do same thing for off-diagonal connectivities
Alpha_offdiag_mean=mean(Alpha_offdiag,2);
S_Alpha_offdiag=(Alpha_offdiag-repmat(Alpha_offdiag_mean,1,size(Alpha_offdiag,2)))*...
    (Alpha_offdiag-repmat(Alpha_offdiag_mean,1,size(Alpha_offdiag,2)))';
Sigma_psi_offdiag=Theta_offdiag*S_Alpha_offdiag*Theta_offdiag';
Sigma_psi_offdiag=.5*(Sigma_psi_offdiag+Sigma_psi_offdiag');
[U D]=eig(Sigma_psi_offdiag);
Theta_offdiag_new=U(:,(K-Q+1):K)*D((K-Q+1):K,(K-Q+1):K)^.5;
Alpha_offdiag_new=regress(squeeze(C{1}(1,2,:)),phi_t'*Theta_offdiag_new);
for i=1:N
    for j1=1:M_i(i)
        for j2=1:M_i(i)
            if j1~=j2
                Alpha_offdiag_new=[Alpha_offdiag_new...
                    regress(squeeze(C{i}(j1,j2,:)),phi_t'*Theta_offdiag_new)];
            end
        end
    end
end
Alpha_offdiag_new=Alpha_offdiag_new(:,2:size(Alpha_offdiag_new,2));
Theta_offdiag=Theta_offdiag_new*cov(Alpha_offdiag_new')^.5;
Alpha_offdiag=cov(Alpha_offdiag_new')^-.5*Alpha_offdiag_new;
Fit_offdiag=phi_t'*Theta_offdiag*Alpha_offdiag;

save smoothing_outputs B Theta_diag C_diag Alpha_diag Fit_diag...
                        C_offdiag Theta_offdiag Alpha_offdiag Fit_offdiag phi_t K Q

% save smoothing coeffients into cell array B for input into MCMC
B=cell(N,1);
ind_diag=0;ind_offdiag=0;
for i=1:N
    B{i}=cell(M_i(i));
    for j1=1:M_i(i)        
        for j2=1:M_i(i)
            if(j1==j2)
                ind_diag=ind_diag+1;
                B{i}{j1,j2}=Alpha_diag(:,ind_diag)';
            end
            
            if(j1~=j2)
                ind_offdiag=ind_offdiag+1;
                B{i}{j1,j2}=Alpha_offdiag(:,ind_offdiag)';
            end
        end
    end
end

%% Save inputs for mcmc algorithm
save real_data S C B N T time M_i subj_keep ind_keep ind_sym Q

