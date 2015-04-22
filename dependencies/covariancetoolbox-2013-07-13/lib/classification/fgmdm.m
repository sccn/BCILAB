function [Ytest,d,C] = fgmdm(COVtest,COVtrain,Ytrain,varargin)

    if isempty(varargin)
        method_mean = 'riemann';
        method_dist = 'riemann';
    else
        method_mean = varargin{1};
        method_dist = varargin{2};
    end
    labels = unique(Ytrain);
    Nclass = length(labels);
    
    % geodesic filtering
    [W,Cg] = fgda(COVtrain,Ytrain,method_mean,{},'shcov',{});
    COVtrain = geodesic_filter(COVtrain,Cg,W(:,1:Nclass-1));
    COVtest = geodesic_filter(COVtest,Cg,W(:,1:Nclass-1));    
    
    [Ytest, d, C] = mdm(COVtest,COVtrain,Ytrain,method_mean,method_dist);