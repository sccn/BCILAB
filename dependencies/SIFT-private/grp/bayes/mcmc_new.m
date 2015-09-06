%function []=mcmc()

clear all; clc;

%% MCMC options
nloop=110; % number of iterations for mcmcm
burnin=10;  % number of initial iterations to discard
thin=10;    % save every nth iteration 
M=4;        %number of group components desired


%% set starting values for Z, S_BAR, SIGMA_S B_BAR SIGMA_B MU
load real_data;
[Z_st S_BAR_st SIGMA_S_st B_BAR_st SIGMA_B_st N_k_st N_k1k2_st]=get_initial_values_new(N,M_i,S,B,Q,M);
Z=Z_st; % Mi x M matrices of group indicators
S_BAR=S_BAR_st; % cluster centroid locations
SIGMA_S=SIGMA_S_st; % cluster centroid variances
B_BAR=B_BAR_st; % group level connectivities
SIGMA_B=SIGMA_B_st; % variances of connectivities
N_k=N_k_st; % number of components that belong to each cluster
N_k1k2=N_k1k2_st; % pairwise counts of group level clusters
MU=ones(M,1)/M; %probabilities of membership in group-level clusters.


%% initialize arrays of posterior draws
keep=round((nloop-burnin)/thin);
Z_array=cell(keep,1);
S_BAR_array=cell(keep,1);
SIGMA_S_array=cell(keep,1);
B_BAR_array=cell(keep,1);
SIGMA_B_array=cell(keep,1);
MU_array=cell(keep,1);
iter_array=0;

%%Hyperparameters
c=10000;                    %'c' in the paper
eta=median(N_k_st);         % 'eta' in paper
SS=eta*mean(SIGMA_S,3);     % 'S' in the paper
p1=1;                       % 'a' in the paper
p2=1;                       % 'b' in the paper


%% run MCMC 

for iter=1:nloop

    % Draw S_BAR (group centroid locations)
    for k=1:M
        mu_s_bar_k=zeros(3,1);
        Sigma_s_bar_k=inv(N_k(k)*inv(SIGMA_S(:,:,k))+inv(SS));
        for i=1:N
            if max(Z{i}(:,k))==1
                j_k=find(Z{i}(:,k)==1);
                for j=1:length(j_k)
                    mu_s_bar_k=mu_s_bar_k+S{i}(j_k(j),:)';
                end
            end
        end
        mu_s_bar_k=Sigma_s_bar_k*inv(SIGMA_S(:,:,k))*mu_s_bar_k;
        S_BAR(k,:)=mvnrnd(mu_s_bar_k,Sigma_s_bar_k);
    end

    
    % Draw SIGMA_S (group centoid covariance matrices)
    for k=1:M
        eta_k=eta+.5*N_k(k);
        SS_k=SS;
        for i=1:N
            if max(Z{i}(:,k))==1
                j_k=find(Z{i}(:,k)==1);
                for j=1:length(j_k)
                    SS_k=SS_k+.5*(S{i}(j_k(j),:)-S_BAR(k,:))'*(S{i}(j_k(j),:)-S_BAR(k,:));
                end
            end
        end
        SS_k=.5*(SS_k+SS_k');
        inv_SS_k=inv(SS_k);
        inv_SS_k=.5*(inv_SS_k+inv_SS_k');
        SIGMA_S(:,:,k)=inv(wishrnd(inv_SS_k,eta_k));
    end  
    
    
    % Draw B_BAR (group mean connectivites)
    for k1=1:M
        for k2=1:M
            mu_b_bar_k1k2=zeros(Q,1);
            Sigma_sq_b_k1k2=1/(1/c+N_k1k2(k1,k2)/SIGMA_B(k1,k2))*eye(Q);
            for i=1:N
                n_k1=sum(Z{i}(:,k1));
                n_k2=sum(Z{i}(:,k2));
                if(and(n_k1>=1,n_k2>=1))
                    j_k1=find(Z{i}(:,k1)==1);
                    j_k2=find(Z{i}(:,k2)==1);
                    for j1=1:n_k1
                        for j2=1:n_k2
                            mu_b_bar_k1k2=mu_b_bar_k1k2+B{i}{j_k1(j1),j_k2(j2)}';
                        end
                    end
                end
            end
            mu_b_bar_k1k2=Sigma_sq_b_k1k2*mu_b_bar_k1k2/SIGMA_B(k1,k2);
            B_BAR{k1,k2}=mvnrnd(mu_b_bar_k1k2,Sigma_sq_b_k1k2);
        end
    end
    
    % Draw SIGMA_B (group connectivity variances)
    for k1=1:M
        for k2=1:M
            p1_k1k2=p1+.5*Q*N_k1k2(k1,k2);
            p2_k1k2=p2;
            for i=1:N
                n_k1=sum(Z{i}(:,k1));
                n_k2=sum(Z{i}(:,k2));
                if(and(n_k1>=1,n_k2>=1))
                    j_k1=find(Z{i}(:,k1)==1);
                    j_k2=find(Z{i}(:,k2)==1);
                    for j1=1:n_k1
                        for j2=1:n_k2
                            p2_k1k2=p2_k1k2+.5*sum((B{i}{j_k1(j1),j_k2(j2)}-B_BAR{k1,k2}).^2);
                        end
                    end
                end
            end
            SIGMA_B(k1,k2)=1/gamrnd(p1_k1k2,1/p2_k1k2);
        end
    end


    % Draw Z (indicators of group membership)   
    ind=randperm(N);
    for i=ind
        ind_i=randperm(M_i(i));
        for j=ind_i
            Log_p_ij=zeros(M,1);
            for k=1:M
                Sigma_s_ijk=SIGMA_S(:,:,k);
                s_ijk=S{i}(j,:)-S_BAR(k,:);
                log_p_s_ijk=-.5*log(det(Sigma_s_ijk))-.5*s_ijk*inv(Sigma_s_ijk)*s_ijk';
                log_p_b_ijk=0;
                for j1=1:M_i(i)
                    if j==j1
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(k,k)))...
                        -.5*sum((B{i}{j,j}-B_BAR{k,k}).^2)/SIGMA_B(k,k);
                    end
                    if and(not(j==j1),find(Z{i}(j1,:)==1)==k)
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(k,k)))...
                        -.5*sum((B{i}{j1,j1}-B_BAR{find(Z{i}(j1,:)==1),find(Z{i}(j1,:)==1)}).^2)/...
                            SIGMA_B(find(Z{i}(j,:)==1),find(Z{i}(j,:)==1));                    
                    end
                    if not(find(Z{i}(j,:)==1)==find(Z{i}(j1,:)==1))
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(k,find(Z{i}(j1,:)==1))))...
                        -.5*sum((B{i}{j,j1}-B_BAR{k,find(Z{i}(j1,:)==1)}).^2)/...
                              SIGMA_B(k,find(Z{i}(j1,:)==1));
                    end
                end
                for j2=1:M_i(i)
                    if not(find(Z{i}(j2,:)==1)==find(Z{i}(j,:)==1))
                        log_p_b_ijk=log_p_b_ijk-.5*Q*log(det(SIGMA_B(find(Z{i}(j2,:)==1),k)))...
                        -.5*sum((B{i}{j2,j}-B_BAR{find(Z{i}(j2,:)==1),k}).^2)/...
                                SIGMA_B(find(Z{i}(j2,:)==1),k);
                    end
                end
                Log_p_ij(k)=MU(k)+log_p_s_ijk+log_p_b_ijk;
            end
            Log_p_ij=Log_p_ij-max(Log_p_ij);
            P_ij=exp(Log_p_ij)/sum(exp(Log_p_ij));
            v=rand;
            Z{i}(j,:)=0*Z{i}(j,:);
            if v<=P_ij(1)
                   Z{i}(j,1)=1;
                   log_p_ij=Log_p_ij(1);
            end
            for k=1:(M-1)
               if(and(v>sum(P_ij(1:k)),v<=sum(P_ij(1:(k+1)))))
                   Z{i}(j,k+1)=1;
                   log_p_ij=Log_p_ij(k);
               end
            end
        end
    end

    % Draw MU (probabilities of clusters)   
    N_k1k2=zeros(M);
    for k1=1:M
        for k2=1:M
            for i=1:N
                n_k1=sum(Z{i}(:,k1));
                n_k2=sum(Z{i}(:,k2));
                if(and(n_k1>=1,n_k2>=1))
                    if not(k1==k2)
                        N_k1k2(k1,k2)=N_k1k2(k1,k2)+n_k1*n_k2;
                    end                
                    if k1==k2
                        N_k1k2(k1,k2)=N_k1k2(k1,k2)+n_k1;
                    end
                end
            end
        end
    end
    N_k=diag(N_k1k2);
    MU=drchrnd(1/M*ones(1,M)+N_k',1);

    
    % Posterior Draw arrays    
    if and(iter>=burnin,thin*round(iter/thin)==iter)
        iter_array=iter_array+1;
        Z_array{iter_array}=Z;
        S_BAR_array{iter_array}=S_BAR;
        SIGMA_S_array{iter_array}=SIGMA_S;
        B_BAR_array{iter_array}=B_BAR;
        SIGMA_B_array{iter_array}=SIGMA_B;
        MU_array{iter_array}=MU;

        %compute and print out intermediate estimates of source locations
        S_BAR_mean=S_BAR_array{1};
        for j=2:iter_array
            S_BAR_mean=S_BAR_mean+S_BAR_array{j};
        end
        S_BAR_mean=S_BAR_mean/iter_array;
        iter
        [S_BAR_st zeros(M,1) S_BAR_mean]
        
        %[sqrt(sum(min(pdist2(S_BAR_true,S_BAR_st)').^2))/M sqrt(sum(min(pdist2(S_BAR_true,S_BAR_mean)').^2))/M]
      
        %compute and print intermediate indicator probs for random subject
        i=discretesample(repmat(1,N,1)/N,1);
        Z_i_mean=Z_array{1}{i};
        for j=2:iter_array
           Z_i_mean=Z_i_mean+Z_array{j}{i};
        end
        Z_st{i}
        Z_i_mean=Z_i_mean/iter_array
    end
        
            
end

save mcmc_outputs nloop M Z_array S_BAR_array SIGMA_S_array B_BAR_array SIGMA_B_array MU_array iter_array

