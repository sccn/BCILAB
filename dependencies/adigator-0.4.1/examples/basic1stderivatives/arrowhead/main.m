% This file uses both MATLAB finite differences as well as adigator in order
% to compute derivatives of the arrowhead function. User can change N to
% see the effects of increasing problem size.
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
N = 100;
numeval = 20;
x = rand(N,1);
% ------------------------------ ADiGator ------------------------------- %
gx = adigatorCreateDerivInput([N, 1],'x'); % Create Deriv Input
gy = adigator('arrowhead',{gx},'arrowhead_dx',adigatorOptions('overwrite',1));
X.f = x;
X.dx = ones(N,1);
tic;
for i = 1:numeval
  y = arrowhead_dx(X);
  dydx = sparse(y.dx_location(:,1),y.dx_location(:,2),y.dx,N,N);
end
adigatortime = toc/numeval;


% -------------------------- Finite Differences ------------------------- %
TOL = 1e-8;% Can change this to make numjac more/less accurate
tic
for i = 1:numeval
  dfdx = numjac(@arrowhead4numjac,0,x,y.f,TOL*ones(N,1),[],0);
end
fdtime = toc/numeval;


% ---------------------- Compressed Finite Differences ------------------ %
S = sparse(y.dx_location(:,1),y.dx_location(:,2),ones(size(y.dx_location(:,1),1),1),N,N);
tic
for i = 1:numeval
  dfdx2 = numjac(@arrowhead4numjac,0,x,y.f,TOL*ones(N,1),[],0,S,[]);
end
cfdtime = toc/numeval;

display(['Average deriv eval time using Finite Differences:           ',num2str(fdtime)]);
display(['Average deriv eval time using ADiGator:                     ',num2str(adigatortime)]);
display(['Average deriv eval time using Compressed Finite Differences:',num2str(cfdtime)]);