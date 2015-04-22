% Generate a set of covariance matrices according to a Wishart distrib.
% Inputs :
%   N : the size of the matrices
%   I : the number of matrices
%   Df : the degee of freedom
%   Sig : the parameter of the wishart disrib 
%
% Outputs :
%   COV : a set of I covariance matrices
%   Sig : the parameter of the wishart disrib 
function [COV Sig] = generate_wishart_set(N,I,Df,Sig)

    if nargin<4
        [Q,~] = qr(randn(N));
        Sig = Q'*diag(5*rand(N,1))*Q;
    end

    COV = zeros(N,N,I);
    for i=1:I
        COV(:,:,i) = wishrnd(Sig,Df);
    end