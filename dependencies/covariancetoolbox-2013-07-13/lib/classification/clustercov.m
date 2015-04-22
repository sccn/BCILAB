function [Ytest d C] = clustercov(Ctest,Ctrain,Nclass,varargin)

    if isempty(varargin)
        method_mean = 'riemann';
        method_dist = 'riemann';
        method_clust = 'complete';
    else
        method_mean = varargin{1};
        method_dist = varargin{2};
        if length(varargin)<=2
            method_clust = 'complete';
        else
            method_clust = varargin{3};
        end
    end
    

    % initialisation
    Ntrial = size(Ctrain,3);

    d = zeros(Ntrial);
    for i=1:Ntrial
        for j=i+1:Ntrial
            d(i,j) = distance(Ctrain(:,:,i),Ctrain(:,:,j),method_dist);
        end
    end
    Y = squareform(d','tovector');
    Z = linkage(Y,method_clust);
    Ytrain = cluster(Z,'maxclust',Nclass);
    
    [Ytest d C] = mdm(Ctest,Ctrain,Ytrain,method_mean,method_dist);
   