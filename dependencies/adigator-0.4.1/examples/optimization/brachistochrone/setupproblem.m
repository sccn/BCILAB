function [probinfo,upperbound,lowerbound,guess, tau] = setupproblem(numintervals,varargin)
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
n = 3;
m = 1;

K = numintervals;
N = numintervals*2+1;
probinfo.numintervals = numintervals;
probinfo.n = n;
probinfo.m = m;

probinfo.k  = 1:2:N-1;
probinfo.kbp1 = probinfo.k+1;
probinfo.kp1    = probinfo.k+2;

probinfo.xind  = reshape(1:n*N,N,n);
probinfo.uind  = reshape(n*N+1:n*N+m*N,N,m);
probinfo.tfind = probinfo.uind(end)+1;


% Set Bounds;
numvars = probinfo.tfind;
upperbound = zeros(numvars,1);
lowerbound = zeros(numvars,1);
% State
upperbound(probinfo.xind) = 10;
lowerbound(probinfo.xind(:,1:2)) = -10;
lowerbound(probinfo.xind(:,3)) = 0;
% Set state boundary conditions
x0ind = probinfo.xind(1,:);
upperbound(x0ind) = 0;
lowerbound(x0ind) = 0;
xfind = probinfo.xind(end,:);
upperbound(xfind(1:2)) = [2;-2];
lowerbound(xfind(1:2)) = [2;-2];

% Control
upperbound(probinfo.uind) = pi;
lowerbound(probinfo.uind) = -pi;
% tf
upperbound(probinfo.tfind) = 5;
lowerbound(probinfo.tfind) = 0;


tau = linspace(0,1,N);
%taum = linspace(0,1,K+1);

% Set Guesses
guess = zeros(numvars,1);
if nargin == 1
  guess(probinfo.xind(:,1)) = linspace(0,2,N);
  guess(probinfo.xind(:,2)) = linspace(0,-2,N);
  guess(probinfo.xind(:,3)) = linspace(0,0,N);
  guess(probinfo.uind) = linspace(-pi,pi,N);
  guess(probinfo.tfind) = 1;
else
  Xold   = varargin{1};
  Uold   = varargin{2};
  tfold  = varargin{3};
  tauold = varargin{4};
  guess(probinfo.xind) = interp1(tauold(:),Xold,tau(:));
  guess(probinfo.uind) = interp1(tauold(:),Uold,tau(:));
  guess(probinfo.tfind) = tfold;
end

% % Build A and B s.t.
% % C = A*X + tf*B*F(X,U)
% 
% % Letting i be first point in interval, j be second, k be third, and let 
% % dt be difference in tau for respective interval
% % C1 = (-1/2*Xi) + (1*Xj) + (-1/2*Xk) + tf*dt*[(-1/8*Fi) +    (0*Fj) +  (1/8*Fk)]
% % C2 =   (-1*Xi) + (0*Xj) +    (1*Xk) + tf*dt*[(-1/6*Fi) + (-4/6*Fj) + (-1/6*Fk)]
% 
% I = repmat(1:K,3,1);
% J = repmat((1:3).',1,K)+repmat(0:2:2*(K-1),3,1);
% A1v = repmat([-1/2; 1; -1/2],1,K);
% A2v = repmat([  -1; 0;    1],1,K);
% 
% deltatau = diag(diff(taum));
% B1v = repmat([-1/8;    0;  1/8],1,K)*deltatau;
% B2v = repmat([-1/6; -4/6; -1/6],1,K)*deltatau;
% 
% A1 = sparse(I,J,A1v,K,N);
% A2 = sparse(I,J,A2v,K,N);
% B1 = sparse(I,J,B1v,K,N);
% B2 = sparse(I,J,B2v,K,N);
% 
% probinfo.A = [A1;A2];
% probinfo.B = [B1;B2];
end