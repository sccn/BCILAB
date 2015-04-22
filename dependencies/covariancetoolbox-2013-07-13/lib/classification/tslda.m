function Ytest = tslda(Ctest,Ctrain,Ytrain,varargin)

    if isempty(varargin)
        method_mean = 'riemann';
        update = 0;
    else
        method_mean = varargin{1};
        update = varargin{2};
    end
    
    labels = unique(Ytrain);
    Nclass = length(labels);
    
    % Tangent space mapping
    C = mean_covariances(Ctrain,method_mean);
    Strain = Tangent_space(Ctrain,C);
    Nelec = size(Strain,1);

    % Regularized LDA
    mu = zeros(Nelec,Nclass);
    Covclass = zeros(Nelec,Nelec,Nclass);

    for i=1:Nclass
        mu(:,i) = mean(Strain(:,Ytrain==labels(i)),2);
        Covclass(:,:,i) = covariances(Strain(:,Ytrain==labels(i)),'shcovft');
    end

    mutot = mean(mu,2);

    Sb = zeros(Nelec,Nelec);    
    for i=1:Nclass
        Sb = Sb+(mu(:,i) - mutot)*(mu(:,i)-mutot)';
    end

    S = mean(Covclass,3);

    [W Lambda] = eig(Sb,S);
    [~, Index] = sort(diag(Lambda),'descend');
    
    W = W(:,Index(1));
    b = W(:,1)'*mutot;

    s = sign(W(:,1)'*mu(:,2)-b);

    % classification
    if update == 1
        C = mean_covariances(Ctest,method_mean);
    end
    Stest = Tangent_space(Ctest,C);
    
    y1 = s*(W(:,1)'*Stest-b);
    Ytest = labels((y1>0) + 1);
    
    
    
   