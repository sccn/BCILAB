% This file uses both MATLAB finite differences as well as adigator in order
% to solve the brusselator ODE using ode15s.
% 
% For this problem, in general, compressed finite differencing seems to be
% the best option as far as time is concerned. However, we compute the
% sparsity matrix using adigator. You can play with N and tspan, but be
% weary that N > 100 will make the non-compressed finite differencing
% extremely slow.

% ----------------------- Problem Set Up -------------------------------- %
fdflag = 1;
N = 100;
tspan = [0; 100];
y0 = [1+sin((2*pi/(N+1))*(1:N)); repmat(3,1,N)];
% Create the Derivative File
tic
gt = adigatorCreateAuxInput([1 1]); % aux input for t
gy = adigatorCreateDerivInput([2*N, 1],'y'); % deriv input for y
gout = adigator('mybrussode',{gt,gy,N},'brussderiv',adigatorOptions('overwrite',1));
adigatorgentime = toc;

% Can build sparsity pattern from output of adigator
S = sparse(gout{1}.deriv.nzlocs(:,1),gout{1}.deriv.nzlocs(:,2),...
  ones(size(gout{1}.deriv.nzlocs,1),1),N*2,N*2);

% ---------------- Solve ODE Using Finite Difference -------------------- %
if fdflag
  tic;
  options = odeset('Vectorized','on'); % Set options
  [t1,y1] = ode15s(@(t,y)mybrussode(t,y,N),tspan,y0,options);   % Solve ODE
  FDtime = toc;
end

% ---------------- Solve ODE Using Compressed Finite Difference --------- %
tic;
options = odeset('Vectorized','on','JPattern',S); % Set options
[t2,y2] = ode15s(@(t,y)mybrussode(t,y,N),tspan,y0,options);     % Solve ODE
CFDtime = toc;

% --------------------- Solve ODE Using ADiGator ------------------------ %
% Now have deriv file generated, but its outputs are not how ode15s likes
% them, so need to write a wrapper file, we put this in
% adigatorodewrapper.m, the code is:
% function [f J] = adigatorodewrapper(t,y)
% Y.f = y;
% Y.dx = ones(size(y));
% dYdt = brussderiv(t,Y);
% f = dYdt.f;
% J = sparse(dYdt.dy_location(:,1),dYdt.dy_location(:,2),...
%   dYdt.dy,dYdt.dy_size(1),dYdt.dy_size(2));
% end
tic
options = odeset('Vectorized','on','Jacobian',@(t,y)adigatorodewrapper(t,y,N),'JPattern',S);
[t,y] = ode15s(@(t,y)mybrussode(t,y,N),tspan,y0,options);     % Solve ODE
adigatortime = toc;
display(['ADiGator deriv file generation time:               ',num2str(adigatorgentime)]);
display(['ODE solve time using ADiGator:                     ',num2str(adigatortime)]);
if fdflag
  display(['ODE solve time using Finite Difference:          ',num2str(FDtime)]);
end
display(['ODE solve time using Compressed Finite Difference: ',num2str(CFDtime)]);