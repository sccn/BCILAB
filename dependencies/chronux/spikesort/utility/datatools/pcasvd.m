function [proj,u,s,v] = pcasvd(data)
%PCASVD            Principal Components Analysis via (mean-subtracted) SVD.
%   PROJ = PCASVD(DATA), where DATA is an M x N matrix, returns the M x N
%   matrix PROJ where the (M,N)th entry is the projection of the Mth row
%   of DATA onto the Nth eigenvector of the covariance matrix COVDATA
%   formed from the rows of DATA (i.e., COVDATA = DATA' * DATA).
%
%   [PROJ,U,S,V] = PCASVD(DATA) also returns matrices U, S, V such that
%   DATA = U * S * V' and PROJ = DATA * V.
%
%   PCASVD always first subtracts the mean from each column of DATA.

% Very simple code -- basically just a macro.
data = removeoffset(data,1,'mean'); % remove mean row
[u,s,v] = svd(data, 0);             % SVD the data matrix
proj = data * v;                    % compute (mean-subtracted) pca projections?