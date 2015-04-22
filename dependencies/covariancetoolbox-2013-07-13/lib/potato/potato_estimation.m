% estimate the center of the ROI and the theeshold given a set of
% covariances matrices.
% [C th dist] = potato_estimation(COV,method_mean,method_th,method_distance,args)
%   Inputs : 
%       * COV : a set of covariance matrices
%       * method_mean : method for center of ROI estimation. 'riemann',
%       'logeuclid','euclidean','geodesic' ... See also MEAN_COVARIANCES
%   See also  MEAN_COVARIANCES, DISTANCE, THREESHOLD_ESTIMATION.

function [C th dist] = potato_estimation(COV,method_mean,method_th,method_distance,args)
    if nargin <2
        method_mean = 'riemann';
        method_th = 'median';
        method_distance = 'riemann';
    end
    if nargin <4
        method_distance = 'riemann';
        args = {};
    end
    if nargin <5
        args = {};
    end
    
    I = size(COV,3);
    
    % estimation of the center of ROI
    C = mean_covariances(COV,method_mean,args);
    
    % distances
    dist = zeros(I,1);
    
    for i=1:I
        dist(i) = distance(COV(:,:,i),C,method_distance);
    end
    
    %threeshold estimation
    th =  threeshold_estimation(dist,method_th);
    
    