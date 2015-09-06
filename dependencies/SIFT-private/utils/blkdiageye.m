function C = blkdiageye(X,k)
% construct block-diagonal matrix with k copies of X on diagonal
% Equivalent to C = kron(eye(k),X) but much faster
%
% Author: Tim Mullen, 2011 (C) SCCN/INC/UCSD

ss = repmat({X},1,k);
C = blkdiag(ss{:});

% alternate form
% ss = repmat('X,',1,k); ss(end)=[];
% eval(sprintf('C = blkdiag(%s);',ss));
