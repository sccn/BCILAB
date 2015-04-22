function th = threeshold_estimation(distance,method)
    if nargin < 2
        method = 'mean';
    end
    
    switch method
        case 'mean'
            th = mean(distance)+2.5*std(distance);
        case 'median'
            th = median(distance)+2.5*1.4826*mad(distance,1);
    end
