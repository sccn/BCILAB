function D = utl_kernelize(X,V,type,scaling,lam)
% Transform some data into a kernel space, given some basis vectors.
% Out-Data = utl_kernelize(In-Data,Basis,Kernel-Type,Scaling,Lambda)
%
% In:
%   In-Data     : input vectors, [NxD], N=#vectors, D=#features
%   Basis       : basis vectors in the same space as the input data, [MxD], M=#basis vectors, D=#features
%   Kernel-Type : type of kernel to apply:
%                  'linear': no kernel
%                  'rbf': gaussian/radial basis function kernel
%                  'poly': polynomial kernel, degree selected by Lambda
%                  'laplace': laplacian kernel
%                  'cauchy': Cauchy kernel
%   Scaling     : scaling of the kernel
%   Lambda      : additional kernel parameter
%
% Out:
%   Out-Data    : data in kernel space
%
% Examples:
%   % map the data into a radial basis function kernel space, using a kernel scale (gamma) of 0.1
%    basis = trials;
%    trials = utl_kernelize(trials,basis,'rbf',0.1);    
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-06

X = X/scaling;
if strcmp(type,'linear')
    D = X;
else
    V = V/scaling;
    dist_squared = @(A,B) repmat(sum(A.^2,2),1,size(B,1)) + repmat(sum(B.^2,2)',size(A,1),1) - 2*(A*B');
    switch type
        case 'rbf'
            D = exp(-dist_squared(X,V));
        case 'poly'
            D = (X*V' + 1).^lam;
        case 'laplace'
            D = exp(-sqrt(dist_squared(X,V)));
        case 'cauchy'
            D = 1./(1+dist_squared(X,V));
        otherwise
            error('unknown kernel type');
    end
end
