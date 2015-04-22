function [Ytest d C] = kmeanscov(Ctest,Ctrain,Nclass,varargin)

    if isempty(varargin)
        method_mean = 'riemann';
        method_dist = 'riemann';
    else
        method_mean = varargin{1};
        method_dist = varargin{2};
    end
    

    % initialisation
    Ntrial = size(Ctrain,3);

    Y_old = zeros(Ntrial,1);
    ix = randperm(Ntrial);
    C = cell(Nclass,1);
    for i=1:Nclass
        C{i} = Ctrain(:,:,ix(i));
    end
    
    d = zeros(Ntrial,Nclass);
    for j=1:Ntrial
        for i=1:Nclass
            d(j,i) = distance(Ctrain(:,:,j),C{i},method_dist);
        end
    end
    
    [~,Y] = min(d,[],2);
    
    % iteration
    while sum(Y~=Y_old)>(Ntrial*0.01)
        Y_old = Y;
        [Y,~,C] = mdm(Ctrain,Ctrain,Y,method_mean,method_dist);        
    end
    
    % classification of test data
    Ntesttrial = size(Ctest,3);
    d = zeros(Ntesttrial,Nclass);
    for j=1:Ntesttrial
        for i=1:Nclass
            d(j,i) = distance(Ctest(:,:,j),C{i},method_dist);
        end
    end
    
    [~,Ytest] = min(d,[],2);

   