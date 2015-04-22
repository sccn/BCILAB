function [COV P1] = covariances_p300(X,Y,method_cov,arg_cov)

    if nargin < 3
        method_cov = 'scm';
        arg_cov = {};
    end
    
    if size(Y)~=size(X(:,:,1))
        P1 = mean(X(:,:,Y==1),3);
    else
        P1 = Y;
    end
    
    X2 = zeros(size(X,1)+size(P1,1),size(X,2),size(X,3));
    for i=1:size(X2,3)
        X2(:,:,i) = cat(1,P1,X(:,:,i));
    end

    % Covariance set
    COV = covariances(X2,method_cov,arg_cov);