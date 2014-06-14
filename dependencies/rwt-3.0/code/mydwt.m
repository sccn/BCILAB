function X = mydwt(X,h,dim,varargin)
%    Y = mydwt(X,h,L,dim);
%
%    Function computes the 1D discrete wavelet transform T for a input signal X using the 
%    scaling filter h.
%
%    Input:
%	    X : signal
%       h : scaling filter
%       Dim : dimension along which to transform (default: first non-singleton dimension)
%       L : number of levels. In the case of a 1D signal, length(x) must be
%           divisible by 2^L; in the case of a 2D signal, the row and the
%           column dimension must be divisible by 2^L. If no argument is
%           specified, a full DWT is returned for maximal possible L.
%
%    Output:
%       Y : the wavelet transform of the signal 

siz = size(X);
if isvector(X)
    X = mdwt(X,h,varargin{:});
else
    % the matrix/tensor case is non-trivial
    if nargin < 3
        dim = find(X~=1,1); end
    % permute X such that dim becomes the first dimension
    if dim ~= 1        
        order = [dim 2:(dim-1) 1 (dim+1):length(siz)];
        X = permute(X,order);
    else
        order = 1:length(siz);
    end
    % matricize X
    X = X(:,:);
    % process each column separately
    for k=1:size(X,2)
        X(:,k) = mdwt(X(:,k),h,varargin{:}); end
    % reverse the previous permutation and matricization
    if dim ~= 1
        X = ipermute(reshape(X,siz(order)),order); end
end
