function X = kronmult(Q, X)
% KRONMULT Efficient Kronecker Multiplication 
%
% Y=kronmult(Q,X) computes 
%     Y = (Q{1} kron Q{2} kron ... kron Q{m})*X
% without ever forming the entire Kronecker product matrix
%     Q{1} kron Q{2} kron ... kron Q{m}.
% This code uses the algorithm from page 394 of Fernandes, et al. 1998,
% JACM 45(3): 381--414 (doi:10.1145/278298.278303).
% The input Q must be a cell array where each Q{i} is a square matrix
% and the input X can be either a vector or a matrix with
% prod(cellfun(@(x) size(x,1),Q)) rows.
%
% See also KRON
%
% Example:
%   ns = [5,4,8,2]; Q=arrayfun(@randn,ns,'UniformOutput',false);
%   x = randn(prod(cellfun(@(x) size(x,1),Q)),1);             % generate data
%   y = kronmult(Q,x);
%   Qfull = 1; for i=1:length(Q), Qfull=kron(Qfull,Q{i}); end % explicit matrix
%   z = Qfull*x;
%   norm(y-z)                                                 % they are equal
%

% Copyright, Stanford University, 2009
% Paul G. Constantine, David F. Gleich

N = length(Q);    % number of matrices
n = zeros(N,1);
nright = 1;
nleft = 1;
for i=1:N-1
    n(i) = size(Q{i},1);
    nleft = nleft*n(i);
end
n(N) = size(Q{N},1);    

for i=N:-1:1
    base = 0;
    jump = n(i)*nright;
    for k=1:nleft
        for j=1:nright
            index1 = base+j;
            index2 = base+j+nright*(n(i)-1);
            X(index1:nright:index2,:) = Q{i}*X(index1:nright:index2,:);
        end
        base = base+jump;
    end
    nleft = nleft/n(max(i-1,1));
    nright = nright*n(i);
end    
