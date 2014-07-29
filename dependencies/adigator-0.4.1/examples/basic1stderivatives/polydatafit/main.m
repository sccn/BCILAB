% This file uses both MATLAB finite differences as well as adigator in order
% to compute derivatives of the fit.m function. The user can change m and n
% to change to problem size.
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
m = 8;
n = 50;
TOL = 1e-5;

x = floor(rand(n,1)*1000)/1000;
d = floor(rand(n,1)*1000)/1000;
numeval = 25;

tic
gx = adigatorCreateDerivInput([n,1],'x');
adigator('fit',{gx,d,m},'fit_x',adigatorOptions('overwrite',1));
gentime = toc;
X.f = x;
X.dx = ones(n,1);
tic
for i = 1:numeval
  p = fit_x(X,d,m);
  % Resulting Jacobian is Full, can just reshape
  dpdx = reshape(p.dx,p.dx_size);
end
adigatortime = toc/numeval;



tic
for i = 1:numeval
  dpdx2 = numjac(@(t,x)fit4numjac(t,x,d,m),0,x,p.f,TOL*ones(n,1),[],0);
end
fdtime = toc/numeval;


fprintf('Derivatives of fit function:\n');
fprintf(['m = %1.0f, n = %1.0f, TOL = ',num2str(TOL),'\n'],m,n);
fprintf(['ADiGator File Generation Time: ',num2str(gentime),'\n']);
fprintf(['ADiGator Average Eval Time:    ',num2str(adigatortime),'\n']);
fprintf(['F Diff Average Eval Time:      ',num2str(fdtime),'\n']);