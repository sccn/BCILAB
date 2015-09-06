function [max_clmass clmembers] = stat_getMaxClusterMass(varargin)


g = arg_define([0 Inf], varargin,...
        arg_norep({'Fd','Observed','Connectivity'},mandatory,[],'Connectivity matrix. This typically contains F-values indicating significant connectivity for each edge. Dimensions must be either [num_nodes x num_nodes x d1 x d2 x ... x dN] where fk is a graph dimension (i.e. time, frequency, etc). Alternately graph can be a [num_features x num_nodes^2] square matrix'), ...
        arg_norep({'Ft','Threshold'},mandatory,[],'Thresholds. This is applied to the connectivity matrix to produce a binary adjacency matrix (connectivity values larger than threshold are kept as edges in the graph). Can be a scalar (constant threshold), a matrix of same dimensionality as Conn containing elementwise thresholds, or a logical matrix of same dimensionality as Conn containing a mask'), ...
        arg({'cdims','ContinuousDims'},[3 4], [],'Continuous dimensions. Indices of dimensions which should be treated as continuous. These dimensions will be clustered using bwlabel.'), ...
        arg({'smetric','SimilarityMetric'},'hamming',{'hamming','jaccard','correlation','euclidean'},'Similarity metric. This is used to compute similarity between the dimensions which are not continuous.'), ...
        arg_subswitch({'clmethod','ClusterMethod'},{'AffinityPropagation'}, ...
        {'AffinityPropagation', ...
            {arg({'prefs','Preferences'},[],[],'Preferences for AP. This is an Nx1 matrix of real numbers called preferences. p(i) indicates the preference that data point i be chosen as a cluster center. A good choice is to set all preference values to the median of the similarity values (default). The number of identified clusters can be increased or decreased  by changing this value accordingly. If p is a scalar, apcluster assumes all preferences are equal to p'), ...
             arg({'options','Options'},[],[],'Name-value pairs containing options for apcluster') ...
            } ...
        }, 'Cluster method to use'), ...
        arg({'minclsize','MinClusterSize'},2,[2 Inf],'Minimum cluster size. See limo_getclustersum for details. Typical values are 2 or 4') ...
        );
%     
% arg_norep({'Fnull','NullDistrib'},[],[],'Null distribution for F. First N-1 dimensions are same as Observed matrix. Last dimension corresponds to samples from distribution (e.g. from bootstrap iterations)'), ...
% arg_norep({'Fnull_th','NullThreshold'},[],[],'Thresholds for null distribution. First N-1 dimensions are same as Threshold matrix. Last dimension corresponds to samples from distribution (e.g. from bootstrap iterations)'), ...
%         
    
% get matrix dims
sz = size(g.Fd);
nd = length(sz);
% continuous dims 
% (to be clustered by bwlabel)
cdims = g.cdims;
% discontinuous dims 
% (to be clustered using apcluster 
% with similarity computed over cdims)
dcdims = setdiff_bc(1:nd,cdims); 

if ~isscalar(g.Ft) && ~isequal(size(g.Ft),sz)
    error('Size of threshold matrix must equal that of Connectivity matrix');
end

% expand scalar threshold to full matrix
if isscalar(g.Ft)
    g.Ft = g.Ft(ones(size(g.Fd))); end
    
% apply thresholding obtain binary mask for Fd 
if islogical(g.Ft)
    mask = g.Ft;
else
    mask = g.Fd > g.Ft;
end

if isempty(dcdims)
   % we will skip clustering of discontinuous dimensions
   % 
   g.Fd = hlp_insertSingletonDim(g.Fd,1);
   mask = hlp_insertSingletonDim(mask,1);
   dcdims = 1;
   cdims = cdims+1;
   
   % reshape matrices to 2D
    if nd>2
        g.Fd = reshape2D(g.Fd,cdims,dcdims);
        if ndims(g.Ft)>1
            mask = reshape2D(mask,cdims,dcdims);
        end
    end
   adj = [];
else
    % Determine neighbors along discontinuous dimensions
    % If there are K discontinuous dimensions of size d1...dK then
    % this procedure forms a [(d1*d2*...*dK) x (d1*d2*...*dK)] adjacency matrix 
    % where each K-tuple has a 1 if the tuples are members of the same cluster with
    % cluster membership determined by the jaccard or other similarity between the vectorized 
    % version of the remaining (continuous) dimensions.
    % reshape matrices to 2D
    if nd>2
        g.Fd = reshape2D(g.Fd,cdims,dcdims);
        if ndims(mask)>1
            mask = reshape2D(mask,cdims,dcdims);
        end
    end

    if ismember_bc(g.smetric,{'jaccard','hamming'})
        % compute similarity matrix on mask
        S = squareform(1-pdist(mask,g.smetric));
        % NaNs occur where the two patterns being compared are zero vectors
        % we consider these to be maximally similar
        S(isnan(S)) = 1; 
        N = size(S,1);
        S(1:N+1:N^2) = 1;
    else
        % compute similarity matrix on original data
        S = squareform(1-pdist(g.Fd,g.smetric));
    end

    clidx = NaN;
    if any(S(:))
        % cluster based on S
        switch g.clmethod.arg_selection
            case 'AffinityPropagation'
                if isempty(g.clmethod.prefs)
                    g.clmethod.prefs = mean(S(:));
                end
                if isempty(g.clmethod.options)
                    g.clmethod.options = {}; end
                [clidx]=apcluster(S,g.clmethod.prefs,g.clmethod.options{:},'maxits',5000);
        end
    end

    v = unique_bc(clidx);
    if all(isnan(clidx)) || length(v) == 1  || all(isnan(v)) == 1
        % there is no special cluster 
        adj = [];
    else
        % remove the most frequent (ie background cluster)  
        % FIXME: mode() perhaps better due to fact that bgcluster may have
        % large similarity in the case where mask is mostly zeros and using
        % hamming distance
%         v = v(v~=mode(clidx)); 
        clmembers = cell(1,length(v));
        clsim     = zeros(1,length(v));
        for s=1:length(v)
            clmembers{s}=find(clidx==v(s));
            % get sum of similarities 
            clsim(s) = sum(abs(S(v(s),clmembers{s})));
        end
        % delete the cluster with minimum sum of similaries
        [~,bgcl] = min(clsim);
        clmembers(bgcl) = [];

        % create QxQ neighbors adjacency matrix where Q = d1 x d2 x ... dn and dk are the sizes of 
        % discontinuous dimensions.
        % adj(i,j) = 1 IFF i and j are in the same cluster, otherwise 0
        adj = false(prod(sz(dcdims)));
        for c=1:size(clmembers,2)
            perms = combnk(clmembers{c},2);
            C = [perms; fliplr(perms)];
            adj(C(:,1),C(:,2)) = true;
        end
        
        % zero out the diagonals
        N = size(adj,1);
        adj(1:N+1:N^2)=0;
    end
    
end

% max cluster mass across all clusters with 'significant' values (i.e. exceeding threshold)
max_clmass = limo_getclustersum(g.Fd,1-mask,adj,g.minclsize,0.5);


end
    
% helper functions
% reshape ND matrix [dcdims(1) x dcdims(2) x ... x cdims(1) x cdims(2) x ... ] to
% 2D matrix [prod(dcdims) x prod(cdims)]
function x = reshape2D(x,cdims,dcdims)
    sz = size(x);
    x = permute(x,[cdims dcdims]);
    x = reshape(x,[prod(sz(cdims)), prod(sz(dcdims))]);
    x = x';
end

