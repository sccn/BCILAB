function Lzz = laghesswrap(z,lambda,probinfo)
% This function builds the Hessian of the Langrangian, where the Lagrangian
% is defined by L = tf + lambda.'*[AA*X + tf*BB*F(X,U)], where AA and BB
% are matrices of size (N-1)*n by N*n, consisting of the matrices A and B,
% respectively, lying on the block diagonals. This function calls second
% derivative files of the dynamics file which computes F(X,U). It then uses
% linear algebra to build the second derivative of L wrt z, Lzz, using the
% outputs of the derivative, Fy, and Fyy, where y = [X U], z = [y tf].
% This is done as follows:
% L     = tf + lambda.'*[AA*X + tf*BB*F(y)]    size         1 by 1
% Ly    = lambda.'*[AA*[1 0] + tf*BB*Fy]       size         1 by N*(n+m)
% Ltf   = 1 + lambda.'*BB*F(y)                 size         1 by 1
% Lyy   = tf*lambda.'*BB*Fyy                   size   N*(n+m) by N*(n+m)
% Ltfy  = lambda.'*BB*Fy                       size         1 by N*(n+m)
% Lytf  = Ltfy.'                               size   N*(n+m) by 1
% Ltftf = 0                                    size         1 by 1
% Lzz   = [Lyy   Lytf]                         size N*(n+m)+1 by N*(n+m)+1
%         [Ltfy     0]
% Where the derivative matrix multiplications must be achieved by unrolling
% dimensions of the derivative matrices.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
n = probinfo.n; % Number of states
m = probinfo.m; % Number of controls
N = probinfo.N; % Number of collocation points

lambda = sparse(lambda.eqnonlin); % Lambda is size 1 by (N-1)*n


X.f  = z(probinfo.xind);  % X  is size N by n
U.f  = z(probinfo.uind);  % U  is size N by m
tf   = z(probinfo.tfind); % tf is size 1 by 1


BB   = probinfo.BB; % BB is of size (N-1)*n by N*n

% Fy is of dimension N by n by N*(n+m), we will need this to compute
% Lytf, so we roll the first two dimensions such that Fy is of size 
% N*n by N*(n+m)
mFy  = N*n;     
nFy  = N*(n+m);
% Fyy is of dimension N by n by N*(n+m) by N*(n+m), we will need this to
% compute Lyy, so we roll the first two dimensions and the last two
% dimensions such that Fyy is of dimension N*n by N*(n+m)*N*(n+m)
mFyy = N*n;
nFyy = N*(n+m)*N*(n+m);
if probinfo.vectders
  % Vectorized Derivatives
  X.dY = ones(size(X.f));
  U.dY = ones(size(U.f));
  
  F    = dynamics_yyvect(X,U,probinfo);
  
  % Build Fy
  iFy  = probinfo.iFy; % previously calculated to roll first 2 dims
  jFy  = probinfo.jFy; 
  Fy   = sparse(iFy,jFy,F.dY,mFy,nFy);

  % Build Fyy
  iFyy = probinfo.iFyy; % previously calculated to roll first 2 dims
  jFyy = probinfo.jFyy; % previously calculated to roll last 2 dims
  Fyy  = sparse(iFyy,jFyy,F.dYdY,mFyy,nFyy);
else
  % Non-Vectorized Derivatives
  X.dy = ones(numel(X.f),1);
  U.dy = ones(numel(U.f),1);
  
  F    = dynamics_yy(X,U,probinfo);
  
  % Build Fy
  % Unroll first two dimensions
  iFy  = sub2ind(F.dy_size(1:2),F.dy_location(:,1),F.dy_location(:,2));
  jFy  = F.dy_location(:,3);
  Fy   = sparse(iFy,jFy,F.dy,mFy,nFy);
  
  % Build Fyy
  % Unroll first two dimensions
  iFyy = sub2ind(F.dydy_size(1:2),F.dydy_location(:,1),F.dydy_location(:,2));
  % Unroll last two dimensions
  jFyy = sub2ind(F.dydy_size(3:4),F.dydy_location(:,3),F.dydy_location(:,4));
  Fyy  = sparse(iFyy,jFyy,F.dydy,mFyy,nFyy);
end

% Lyy
Lyy = tf*lambda.'*BB*Fyy;             % 1 by N*(n+m)*N*(n+m)
Lyy = reshape(Lyy,[N*(n+m) N*(n+m)]); % N*(n+m) by N*(n+m)

% Ltfy
Ltfy = lambda.'*BB*Fy;  % 1 by N*(n+m)
Lytf = Ltfy.';          % N*(n+m) by 1

% Lzz
Lzz = [Lyy Lytf;Ltfy 0];