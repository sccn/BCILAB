function MCMC_State = grp_hbmm_init(varargin)
% Initialize MCMC Sampler
% Inputs: 
% 
% Outputs:
%
% MCMC_State
%     .Z         : Mi x M matrices of group indicators
%     .S_BAR     : cluster centroid locations
%     .SIGMA_S   : cluster centroid variances
%     .B_BAR     : group level connectivities
%     .SIGMA_B   : variances of connectivities
%     .N_k       : number of components that belong to each cluster
%     .N_k1k2    : pairwise counts of group level clusters
%
% Author: Tim Mullen and Wes Thompson, SCCN/INC/UCSD, 2010-2014

arg_define([0 2], varargin, ...
    arg_norep({'B','Connectivity'},mandatory,[],'Connectivity matrices. B{i} is an (M_i x M_i x Q) matrix of (possibly time-varying) connectivity values or basis coefficients for the ith subject. M_i is the number of components for the ith subject.'), ...
    arg_norep({'S','DipoleLocations'},mandatory,[],'Dipole locations. S{i} is an (M_i x 3) matrix of [X Y Z] dipole locations for the ith subject. M_i is the number of components for the ith subject.'), ...  
    arg_subswitch({'initVals','InitialValues'},'clustering', ...
        {'clustering', { ...
            arg_subswitch({'clmethod','ClusterMethod'},'KMeans', ...
                {'AffinityPropagation', ...
                    {...
                    arg({'distance','DistanceMetric'},'euclidean',{'euclidean','seuclidean','cityblock','minkowski','chebychev','mahalanobis','cosine','correlation','spearman','hamming','jaccard'},'Distance metric. See ''doc pdist'' for details.'), ...
                    arg({'prefs','Preferences'},[],[],'Preferences for AP. Ignored if NumClusters is specified. This is an [num_observations x 1] vector of real numbers called preferences. p(i) indicates the preference that data point i be chosen as a cluster center. A good choice is to set all preference values to the median of the similarity values (default). The number of identified clusters can be increased or decreased  by changing this value accordingly. If p is a scalar, apcluster assumes all preferences are equal to p'), ...
                    arg({'numclust','NumClusters'},[],[],'Number of clusters desired. Leave empty to determine automatically (based on Preferences)'), ...
                    arg({'options','Options'},[],[],'Name-value pairs containing additional options for apcluster. Type ''doc apcluster'' for details.') ...
                    } ...
                'KMeans', ...
                    {...
                    arg({'distance','DistanceMetric'},'sqEuclidean',{'Hamming','correlation','cityblock','sqEuclidean','cosine'},'Distance metric. See ''doc kmeans'' for details.'), ...
                    arg({'numclust','NumClusters'},[],[],'Number of clusters desired.'), ...
                    arg({'options','Options'},[],[],'Name-value pairs containing additional options for kmeans. Type ''doc kmeans'' for details.') ...
                    } ...
                'GMM', ...
                    {...
                    arg({'numclust','NumClusters'},[],[],'Number of clusters desired.'), ...
                    arg({'options','Options'},[],[],'Name-value pairs containing additional options for gmmfit. Type ''doc gmdistribution.fit'' for details.') ...
                    } ...
                },'Cluster method to use'), ...
            }, ...
        },'MCMC initial clustering (group-level means and covariances). This can be determined automatically via clustering or supplied by the user.') ...
    );

if ~iscell(B)
    error('B must be a cell array'); end
if ~iscell(S)
    error('S must be a cell array'); end
if length(B)~=length(S)
    error('B and S must be the same length (same number of subjects)'); end

% number of subjects
N = length(B);

% --------------------------
% perform clustering
if strcmp(initVals.arg_selection,'clustering')
    opts = initVals.clmethod.options;
    if isempty(opts) && ~iscell(opts), opts = {}; end
    % concatenate all dipole observations into [NN x 3] matrix
    S_clust = cell2mat(S);
    switch initVals.clmethod.arg_selection
        case 'AffinityPropagation'
            % Find clusters using Affinity Propagation
            % ------------------------------------------
            dmetric = initVals.clmethod.distance;
            if ismember_bc(dmetric,{'jaccard','hamming'})
                % compute similarity matrix on mask
                SM = squareform(1-pdist(S_clust,dmetric));
                % NaNs occur where the two patterns being compared are zero vectors
                % we consider these to be maximally similar
                SM(isnan(SM)) = 1; 
                H = size(SM,1);
                SM(1:H+1:H^2) = 1;
            else
                % compute similarity matrix on original data
                S = squareform(1-pdist(S_clust,dmetric));
            end
            if isempty(opts)
                opts = {'maxits',1000}; end
            if ~isempty(initVals.clmethod.numclust)
                % find approximately numclust clusters
                [idx]=apclusterK(S,initVals.clmethod.numclust);
            else
                % find data-determined number of clusters
                if isempty(initVals.clmethod.prefs)
                    initVals.clmethod.prefs = median(S);
                end
                [idx]=apcluster(S,initVals.clmethod.prefs,opts{:});
            end
            uidx = unique_bc(idx);
            % compute cluster centroids
            cc = arrayfun(@(ci)mean(S(idx==ci,:),1), uidx, 'UniformOutput',false);
            cc = cell2mat(cc);
            % recode idx as [1 2 ... M] where M is # clusters
            idxb = idx;
            for i=1:length(uidx)
                idx(idxb==uidx(i)) = i; end
        case 'KMeans'
            % Find clusters using K-Means
            % ------------------------------------------
            [idx cc]=kmeans(S_clust,initVals.clmethod.numclust,     ...
                            'distance',initVals.clmethod.distance,  ...
                            opts{:});
        case 'GMM'
            % Find clusters using Gaussian Mixture Model
            % ------------------------------------------
            gmfit=gmdistribution.fit(S_clust,initVals.clmethod.numclust,opts{:});
            cc=gmfit.mu;
            idx=cluster(gmfit,S_clust);
    end
end

% number of desired clusters
M = length(unique_bc(idx));
% number of dipoles for each subject
M_i = cellfun(@(B_i) size(B_i,1),B);
% length of connectivity sequence (same for all subjects)
Q = size(B{1},3);

% --------------------------
% re-order dipoles to match clusters across subjects
cc_sorted=cc;
xx=sort(cc(:,1));
ord=zeros(M,1);
for k=1:M
    cc_sorted(k,:)=cc(cc(:,1)==xx(k),:);
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
% FIXME: replace with this:
% SIGMA_S = arrayfun(@(ci)cov(S_clust(idx==ci,:),1), idx, 'UniformOutput',false);
SIGMA_S=zeros(3,3,M);
for k=1:M
    SIGMA_S(:,:,k)=cov(S_clust(idx==k,:));
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
B_BAR =zeros(M,M,Q);
for k1=1:M
    for k2=1:M
        B_BAR(k1,k2,:)=zeros(1,Q);
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
                        B_BAR(k1,k2,:)=squish(B_BAR(k1,k2,:))+squish(B{i}(ind1(j1),ind2(j2),:));
                    end
                end
            end
        end
        B_BAR(k1,k2,:)=B_BAR(k1,k2,:)/N_k1k2(k1,k2);
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
                        Sigma_sq_b_k1k2=Sigma_sq_b_k1k2+sum((B{i}(ind1(j1),ind2(j2),:)-...
                            B_BAR(k1,k2,:)).^2)/Q;
                    end
                end
            end
        end
        SIGMA_B(k1,k2)=Sigma_sq_b_k1k2/(1+N_k1k2(k1,k2));
    end
end

% --------------------------
% build output structure
%
% MCMC_State
%     .Z         : Mi x M matrices of group indicators
%     .S_BAR     : cluster centroid locations
%     .SIGMA_S   : cluster centroid variances
%     .B_BAR     : group level connectivities
%     .SIGMA_B   : variances of connectivities
%     .N_k       : number of components that belong to each cluster
%     .N_k1k2    : pairwise counts of group level clusters
varnames = {'Z','S_BAR','SIGMA_S','B_BAR','SIGMA_B','N_k','N_k1k2'};
for i = 1:length(varnames)
    MCMC_State.(varnames{i}) = eval(varnames{i}); 
end
MCMC_State.initstate = true;
