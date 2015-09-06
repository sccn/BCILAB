function X = hlp_insertSingletonDim(X,dim)
% insert a singleton dimension into dimension 'dim' of matrix X
% author: Tim Mullen, 2011

    nd = ndims(X);
    X = permute(X,unique_bc([1:dim-1 nd+1 dim:nd]));
    