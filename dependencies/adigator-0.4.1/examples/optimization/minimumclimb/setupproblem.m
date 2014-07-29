function [probinfo,upperbound,lowerbound,guess, tau] = setupproblem(numintervals,varargin)
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
n = 3; %number of states
m = 1; %number of controls
K = numintervals;
N = K*2+1;
probinfo.numintervals = K;
probinfo.n = n;
probinfo.m = m;
probinfo.N = N;

probinfo.xind  = reshape(1:n*N,N,n);
probinfo.uind  = reshape(n*N+1:n*N+m*N,N,m);
probinfo.tfind = probinfo.uind(end)+1;


% Some problem information for use in bounds
feettometer = .3048;
hscale = 1000;
tscale = 200;

probinfo.hscale = hscale;
probinfo.tscale = tscale;

h0 = 0;
hf = 20;
v0 = 26;
vf = 59;
fpa0 = 0;
fpaf = 0;

hmin =  0*feettometer       /hscale;
hmax =  69000*feettometer   /hscale;
vmin =  1*feettometer       /hscale*tscale;
vmax =  2000*feettometer    /hscale*tscale;
fpamin = -90/180*pi;
fpamax = -fpamin;
umin = -10;
umax =  10;
tfmin = 100 /tscale;
tfmax = 350 /tscale;

% Set Bounds
numvars = probinfo.tfind;
upperbound = zeros(numvars,1);
lowerbound = zeros(numvars,1);
% State
upperbound(probinfo.xind(:,1)) = hmax; % Altitude
lowerbound(probinfo.xind(:,1)) = hmin;
upperbound(probinfo.xind(:,2)) = vmax; % Speed
lowerbound(probinfo.xind(:,2)) = vmin; 
upperbound(probinfo.xind(:,3)) = fpamax; % flight path angle
lowerbound(probinfo.xind(:,3)) = fpamin;

% Set state boundary conditions
x0ind = probinfo.xind(1,:);   % Initial state
upperbound(x0ind) = [h0; v0; fpa0];
lowerbound(x0ind) = [h0; v0; fpa0];
xfind = probinfo.xind(end,:); % Final State
upperbound(xfind) = [hf; vf; fpaf];
lowerbound(xfind) = [hf; vf; fpaf];

% Control
upperbound(probinfo.uind) = umax;
lowerbound(probinfo.uind) = umin;
% tf
upperbound(probinfo.tfind) = tfmax;
lowerbound(probinfo.tfind) = tfmin;


tau  = linspace(0,1,N);
taum = linspace(0,1,K+1);

% Set Guesses
guess = zeros(numvars,1);
if nargin == 1
  guess(probinfo.xind(:,1)) = linspace(hmin,hmax,N);
  guess(probinfo.xind(:,2)) = linspace(vmin,vmax,N);
  guess(probinfo.xind(:,3)) = linspace(0,0,N);
  guess(probinfo.uind)      = linspace(0,0,N);
  guess(probinfo.tfind)     = 1;
else
  Xold   = varargin{1};
  Uold   = varargin{2};
  tfold  = varargin{3};
  tauold = varargin{4};
  guess(probinfo.xind)  = interp1(tauold(:),Xold,tau(:),'spline');
  guess(probinfo.uind)  = interp1(tauold(:),Uold,tau(:),'spline');
  guess(probinfo.tfind) = tfold;
end

% Build A and B s.t.
% C = A*X + tf*B*F(X,U)

% Letting i be first point in interval, j be second, k be third, and let 
% dt be difference in tau for respective interval
% C1 = (-1/2*Xi) + (1*Xj) + (-1/2*Xk) + tf*dt*[(-1/8*Fi) +    (0*Fj) +  (1/8*Fk)]
% C2 =   (-1*Xi) + (0*Xj) +    (1*Xk) + tf*dt*[(-1/6*Fi) + (-4/6*Fj) + (-1/6*Fk)]

I = repmat(1:K,3,1);
J = repmat((1:3).',1,K)+repmat(0:2:2*(K-1),3,1);
A1v = repmat([-1/2; 1; -1/2],1,K);
A2v = repmat([  -1; 0;    1],1,K);

deltatau = diag(diff(taum));
B1v = repmat([-1/8;    0;  1/8],1,K)*deltatau;
B2v = repmat([-1/6; -4/6; -1/6],1,K)*deltatau;

Ia  = [I(:);I(:)+K];
Ja  = [J(:);J(:)];
Ava = [A1v(:);A2v(:)];
Bva = [B1v(:);B2v(:)];

probinfo.A = sparse(Ia,Ja,Ava,2*K,N);
probinfo.B = sparse(Ia,Ja,Bva,2*K,N);

% Big B is just B repmatted on the diagonal n times
% Rows inds get added by 2K*(0:n-1)
Ib  = repmat(Ia,1,n)+repmat(2*K*(0:n-1),numel(Ia),1);
% Cols inds get added by  N*(0:n-1)
Jb  = repmat(Ja,1,n)+repmat(  N*(0:n-1),numel(Ja),1);
Bvb = repmat(Bva,1,n);

probinfo.BB  = sparse(Ib,Jb,Bvb,2*K*n,N*n);


% Constants for Interpolation in this problem
CoF      = zeros(10,6);
CoF(1,:) = [2.61059846050e-2;
            -8.57043966269e-2;
            1.07863115049e-1;
            -6.44772018636e-2;
            1.64933626507e-2;
            0];
CoF(2,:) = [1.37368651246e0;
            -4.57116286752e0;
            5.72789877344e0;
            -3.25219000620e0;
            7.29821847445e-1;
            0];
CoF(3,:) = [1.23001735612e0;
            -2.97244144190e0;
            2.78009092756e0;
            -1.16227834301e0;
            1.81868987624e-1;
            0];
CoF(4,:) = [1.42392902737e1;
            -3.24759126471e1;
            2.96838643792e1;
            -1.33316812491e1;
            2.87165882405e0;
            -2.27239723756e-1];
CoF(5,:) = [0.11969995703e6;
            -0.14644656421e5;
            -0.45534597613e3;
            0.49544694509e3;
            -0.46253181596e2;
            0.12000480258e1];
CoF(6,:) = [-0.35217318620e6;
            0.51808811078e5;
            0.23143969006e4;
            -0.22482310455e4;
            0.20894683419e3;
            -0.53807416658e1];
CoF(7,:) = [0.60452159152e6;
            -0.95597112936e5;
            -0.38860323817e4;
            0.39771922607e4;
            -0.36835984294e3;
            0.94529288471e1];
CoF(8,:) = [-0.43042985701e6;
            0.83271826575e5;
            0.12357128390e4;
            -0.30734191752e4;
            0.29388870979e3;
            -0.76204728620e1];
CoF(9,:) = [0.13656937908e6;
            -0.32867923740e5;
            0.55572727442e3;
            0.10635494768e4;
            -0.10784916936e3;
            0.28552696781e1];
CoF(10,:) = [-0.16647992124e5;
            0.49102536402e4;
            -0.23591380327e3;
            -0.13626703723e3;
            0.14880019422e2;
            -0.40379767869e0];

probinfo.CONSTANTS.CoF = CoF;

end