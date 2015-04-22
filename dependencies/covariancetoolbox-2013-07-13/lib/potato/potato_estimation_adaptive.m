% estimate the center of the ROI and the theeshold given a set of
% covariances matrices.
function [C th dist art Ci thi] = potato_estimation_adaptive(COV,method_mean,method_distance,window,C)
    I = size(COV,3);
    if nargin <2
        method_mean = 'riemann';
        method_distance = 'riemann';
        window = I;
    end
    if nargin <3
        method_distance = method_mean;
        window = I;
    end
    if nargin < 4
        window = I;
    end
    if nargin <5  
        C = COV(:,:,1);
    end
    % distances
    dist = zeros(I,1);

    % artefact[C th dist art Ci thi]
    art = zeros(I,1);
        
    %mean distance
    mu = zeros(I,1);
    
    %std of distance
    std = ones(I,1);
    
    k=1;    
    for i=2:I
        dist(i) = distance(COV(:,:,i),C,method_distance);
        if (dist(i) < mu(i-1)+2.5*std(i-1))||(i<10);
            k=k+1;
            alpha = min([k window]);
            mu(i) = ((alpha-1)/alpha)*mu(i-1)+(1/alpha)*dist(i);
            std(i) = sqrt(((alpha-1)/alpha)*(std(i-1)^2)+(1/alpha)*(dist(i)-mu(i))^2);
            C = geodesic(C,COV(:,:,i),1/alpha,method_mean);
        else
            mu(i) = mu(i-1);
            std(i) = std(i-1);
            art(i) = 1;
        end
        Ci(:,:,i) = C;
        thi(i) = mu(i)+2.5*std(i);
    end
    
    th = mu(i)+2.5*std(i);
    
    
    