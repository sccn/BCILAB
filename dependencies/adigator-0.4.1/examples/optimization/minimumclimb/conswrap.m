function [C,Ceq,JC,JCeq] = conswrap(z,probinfo)
% This function builds the derivatives of C wrt z by calling the dynamics
% derivative file which computes the derivatives of F wrt y, where %
% y = [X U], z = [y tf]. This is done using the fact that C is written as 
% C = A*X + tf*B*F(X,U). Thus, we only take the derivatives of F and use
% linear algebra to build the derivative of C. Note: C is a matrix of size
% (N-1) by n, so we unroll this to create a vector of constraints, that is,
% C = C(:), and when we say C, we mean Ceq since we only have equality
% constraints. We build this as follows:
% C     = A*X + tf*B*F(X,U)
% dCdX  = A*dXdX + tf*B*dFdX
% dCdU  = tf*B*dFdU
% dCdtf = B*F
% dCdz  = [dCdX dCdU dCdtf]
% Where the derivative matrix multiplications must be achieved by
% unrolling dimensions of the derivative matrices.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
n = probinfo.n; % Number of states
m = probinfo.m; % Number of controls
N = probinfo.N; % Number of collocation points

A = probinfo.A; % A is size (N-1) by N
B = probinfo.B; % B is size (N-1) by N

% C = A*X + B*F(X,U) is size (N-1) by n

C  = [];
JC = [];
if nargout == 2
  % Doesn't want derivatives
  X  = z(probinfo.xind);
  U  = z(probinfo.uind);
  tf = z(probinfo.tfind);
  F   = dynamics(X,U,probinfo);
  Ceq = probinfo.A*X + tf.*probinfo.B*F; 
else
  X.f  = z(probinfo.xind);  % X  is size N by n
  U.f  = z(probinfo.uind);  % U  is size N by m
  tf   = z(probinfo.tfind); % tf is size 1 by 1
  
  % dFdy is of dimension N by n by N*(n+m), in this section we unroll
  % the first two dimensions such that our dFdy is of dimension 
  % N*n by N*(n+m)
  mFy  = N*n;
  nFy  = N*(n+m);
  if probinfo.vectders
    % Vectorized Derivatives
    X.dY = ones(size(X.f));
    U.dY = ones(size(U.f));
    
    F    = dynamics_yvect(X,U,probinfo); % F    is size N by n
    iFy  = probinfo.iFy; % previously calculated to roll 1st 2 dims
    jFy  = probinfo.jFy; % (by adigatorProjectVectLocs)
    dFdy = sparse(iFy,jFy,F.dY,mFy,nFy);
  else
    % Non-Vectorized Derivatives
    X.dy = ones(numel(X.f),1);
    U.dy = ones(numel(U.f),1);
    
    F    = dynamics_y(X,U,probinfo);     % F    is size N by n
    % Roll first two dimensions
    iFy  = sub2ind(F.dy_size(1:2),F.dy_location(:,1),F.dy_location(:,2));
    jFy  = F.dy_location(:,3);
    dFdy = sparse(iFy,jFy,F.dy,mFy,nFy); 
  end
  
  % Extract and reshape dFdX, dFdU, and build dXdX
  % Here we switch from 
  % dXdX in N*n by N*n to N by n*N*n
  % dFdX in N*n by N*n to N by n*N*n
  % dFdU in N*n by N*m to N by n*N*m
  dXdX = reshape(sparse(1:N*n,1:N*n,1,N*n,N*n),[N,n*N*n]);
  dFdX = reshape(dFdy(:,probinfo.xind(:)),     [N,n*N*n]);
  dFdU = reshape(dFdy(:,probinfo.uind(:)),     [N,n*N*m]);
  % dXdX, dFdX and dFdU are now in their reshaped, unrolled form.
  
  
  % The output of this is derivative of C wrt z, where z = [X(:); U(:); tf]
  % It is given in the unrolled form, thus
  % C    is size N-1 by 1
  % dCdz is size N-1 by (m+n)*N+1
  
  % Constraints
  Ceq = A*X.f + tf*B*F.f;             % (N-1)   by n
  Ceq = reshape(Ceq,    [(N-1)*n,1]); % (N-1)*n by 1
  
  % Derivative wrt x
  dCdX  = A*dXdX + tf.*B*dFdX;          % (N-1)   by n*N*n
  dCdX  = reshape(dCdX, [(N-1)*n,n*N]); % (N-1)*n by N*n
  
  % Derivative wrt u
  dCdU  = tf.*B*dFdU;                   % (N-1)   by n*N*m
  dCdU  = reshape(dCdU, [(N-1)*n,m*N]); % (N-1)*n by N*m
  
  % Derivative wrt tf
  dCdtf = B*F.f;                        % (N-1)   by n
  dCdtf = reshape(dCdtf,[(N-1)*n,1]);   % (N-1)*n by 1
  
  % Build dCdz
  dCdz = [dCdX dCdU dCdtf];  % (N-1)*n by N*(m+n) + 1
  
  % fmincon takes the transpose of it.
  JCeq = dCdz.';
end
