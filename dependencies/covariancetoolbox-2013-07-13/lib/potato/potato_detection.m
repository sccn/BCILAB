function [art, d] = potato_detection(COV,C,th,method_dist)

    if nargin <4
        method_dist = 'riemann';
    end

    I = size(COV,3);
    d = zeros(I,1);

    for i=1:I
        d(i) = distance(C,COV(:,:,i),method_dist);
    end

    art = d>th;