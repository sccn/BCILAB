function [Ytest d C] = mdm(COVtest,COVtrain,Ytrain,varargin)
    
    if isempty(varargin)
        method_mean = 'riemann';
        method_dist = 'riemann';
    else
        method_mean = varargin{1};
        method_dist = varargin{2};
    end
    
    labels = unique(Ytrain);
    Nclass = length(labels);
    C = cell(Nclass,1);
    
    % estimation of center
    for i=1:Nclass
        C{i} = mean_covariances(COVtrain(:,:,Ytrain==labels(i)),method_mean);
    end

    % classification
    NTesttrial = size(COVtest,3);
    
    d = zeros(NTesttrial,Nclass);
    for j=1:NTesttrial
        for i=1:Nclass
            d(j,i) = distance(COVtest(:,:,j),C{i},method_dist);
        end
    end
    
    [~,ix] = min(d,[],2);
    Ytest = labels(ix);