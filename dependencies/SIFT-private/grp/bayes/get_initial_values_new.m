function [Z S_BAR SIGMA_S B_BAR SIGMA_B N_k N_k1k2]=get_initial_values_new(N,M_i,S,B,Q,M,k)
% Z         : Mi x M matrices of group indicators
% S_BAR     : cluster centroid locations
% SIGMA_S   : cluster centroid variances
% B_BAR     : group level connectivities
% SIGMA_B   : variances of connectivities
% N_k       : number of components that belong to each cluster
% N_k1k2    : pairwise counts of group level clusters

% initial k-means clustering step
if k==1
    S_clust=S{1};
    for i=2:N
        S_clust=[S_clust; S{i}];
    end
    [idx cc]=kmeans(S_clust,M);
end
% ... or GMM clustering
if not(k==1)
    S_clust=S{1};
    for i=2:N
        S_clust=[S_clust; S{i}];
    end
    gmfit=gmdistribution.fit(S_clust,M);
    cc=gmfit.mu;
    idx=cluster(gmfit,S_clust);
end

% --------------------------
% re-order dipoles to match clusters across subjects
cc_sorted=cc;
xx=sort(cc(:,1));
ord=zeros(M,1);
for k=1:M
    cc_sorted(k,:)=cc(find(cc(:,1)==xx(k)),:);
    ord(k)=find(cc(:,1)==xx(k));
end
idx_sorted=idx;
for m=1:length(idx)
    idx_sorted(m)=find(ord==idx(m));
end
cc=cc_sorted;
idx=idx_sorted;

% --------------------------
% dipole location centroids
S_BAR=cc;

% --------------------------
% compute dipole location covariance matrices
SIGMA_S=zeros(3,3,M);
for k=1:M
    SIGMA_S(:,:,k)=cov(S_clust(find(idx==k),:));
end

% --------------------------
% initialize indicator variables Z
% (Z{i}(j,k) = 1 IFF comp_j from subject i belongs to cluster k, otherwise 0
Z=cell(N,1);
ind=0;
for i=1:N
    Z{i}=zeros(M_i(i),M);
    for j=1:M_i(i)
        ind=ind+1;
        for k=1:M
            if idx(ind)==k
                Z{i}(j,k)=1;
            end
        end
    end
end

% --------------------------
% compute between-cluster connectivity mean
N_k1k2=zeros(M);
B_BAR=cell(M);
for k1=1:M
    for k2=1:M
        B_BAR{k1,k2}=zeros(1,Q);
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
                ind1=find(Z{i}(:,k1)==1);
                ind2=find(Z{i}(:,k2)==1);
                for j1=1:n_k1
                    for j2=1:n_k2
                        B_BAR{k1,k2}=B_BAR{k1,k2}+B{i}{ind1(j1),ind2(j2)};
                    end
                end
            end
        end
        B_BAR{k1,k2}=B_BAR{k1,k2}/N_k1k2(k1,k2);
    end
end
N_k=diag(N_k1k2);

% --------------------------
% compute between-cluster connectivity covariance matrices
SIGMA_B=zeros(M);
for k1=1:M
    for k2=1:M
        Sigma_sq_b_k1k2=1;
        for i=1:N
            n_k1=sum(Z{i}(:,k1));
            n_k2=sum(Z{i}(:,k2));
            if(and(n_k1>=1,n_k2>=1))
                ind1=find(Z{i}(:,k1)==1);
                ind2=find(Z{i}(:,k2)==1);
                for j1=1:n_k1
                    for j2=1:n_k2
                        Sigma_sq_b_k1k2=Sigma_sq_b_k1k2+sum((B{i}{ind1(j1),ind2(j2)}-...
                            B_BAR{k1,k2}).^2)/Q;
                    end
                end
            end
        end
        SIGMA_B(k1,k2)=Sigma_sq_b_k1k2/(1+N_k1k2(k1,k2));
    end
end
