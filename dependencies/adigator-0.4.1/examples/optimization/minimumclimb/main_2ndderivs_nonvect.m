% This function solves the minimum time to climb problem by supplying first
% and second derivatives without advantage of the vectorized nature of the
% problem. A derivative file for the dynamics is created as dynamics_yvect
% and these derivatives are used to build the constraint jacobian within
% the file conswrap.m. Also, a second derivative file for the dynamics is
% created as dynamics_yyvect, these derivatives, along with the first
% derivatives, are used to build the lagrangian hessian within the file
% laghesswrap.m.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
options = optimset('Algorithm','interior-point','MaxFunEvals',50000,...
  'GradObj','on','GradConstr','on','Display','iter','MaxIter',1000);
time = zeros(4,1);
numintervals = [10,20,40];
for i = 1:3
tic;
% ---------------------- Set Up the Problem ----------------------------- %
if i == 1
[probinfo,upperbound,lowerbound,guess,tau] = setupproblem(numintervals(i));
else
% If on second or third mesh, use previous solution as initial guess
[probinfo,upperbound,lowerbound,guess,tau] = ...
  setupproblem(numintervals(i),X,U,fval,tau);
end
probinfo.vectders = 0;
options = optimset(options,'Hessian','user-supplied','HessFcn',...
  @(z,lambda)laghesswrap(z,lambda,probinfo));
% ------------------- Create Dynamics Derivative File ------------------- %
% NOTE: fmincon is given the function conswrap which builds the
% Jacobian by calling dynamics_y

% Inputs to dynamics file are X,U,probinfo
% X is size N by n
% U is size N by m
% Let y = [X(:);U(:)] - take derivatives of dynamics file wrt y
% Note: our decision vector is z = [y;tf], and dynamics file is not a
% function of tf, thus dFdz = [dFdy zeros((N+1)*n,1]
n = probinfo.n; m = probinfo.m; N = probinfo.N; 

gX = adigatorCreateDerivInput([N n],...
  struct('vodname','y','vodsize',[(m+n)*N 1],...
  'nzlocs',[probinfo.xind(:) probinfo.xind(:)]));
gU = adigatorCreateDerivInput([N m],...
  struct('vodname','y','vodsize',[(m+n)*N 1],...
  'nzlocs',[(1:m*N).' probinfo.uind(:)]));
adigator('dynamics',{gX,gU,probinfo},'dynamics_y',...
  adigatorOptions('overwrite',1));
% Make new vars for gX,gU
gX = struct('f',gX,'dy',ones(n*N,1));
gU = struct('f',gU,'dy',ones(m*N,1));
adigator('dynamics_y',{gX,gU,probinfo},'dynamics_yy',...
  adigatorOptions('overwrite',1,'comments',0));
% Also note, we must create new derivative files each time we change
% meshes

% --------------------------- Call fmincon ------------------------------ %
[z,fval] = ...
  fmincon(@(x)basic_objwrap(x,probinfo),guess,[],[],[],[],...
  lowerbound,upperbound,@(x)conswrap(x,probinfo),options);

% --------------------------- Extract Solution -------------------------- %
t = tau(:).*z(end);
X = z(probinfo.xind);
U = z(probinfo.uind);
time(i) = toc;
disp(z(end))
% ---------------------------- Plot Solution ---------------------------- %
figure(i+3);
subplot(2,1,1)
plot(t,X,'-o');
xlabel('time')
ylabel('states')
title(sprintf(['Mesh %1.0f Solution in ',num2str(time(i)),'s'],i))
subplot(2,1,2);
plot(t,U,'-o');
xlabel('time')
ylabel('control')

end

fprintf(['Total Time Supplying Second Derivatives (non-vectorized): ',...
  num2str(sum(time)),'\n']);