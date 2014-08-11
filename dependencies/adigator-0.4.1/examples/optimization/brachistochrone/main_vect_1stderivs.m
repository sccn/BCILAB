% This function solves the brachistochrone problem by supplying first
% derivatives and takes advantage of the vectorized nature of the
% problem.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
options = optimset('Algorithm','interior-point','MaxFunEvals',50000,...
  'GradObj','on','GradConstr','on','Display','iter');
time = zeros(4,1);
numintervals = [5,10,20,40];
tic
% -------------- Vectorized Derivatives of Dynamics File ---------------- %
% We want to take the derivatives of the dynamics file since it is a file
% written in a vectorized manner which computes F(Z(t)), where
% Y(t)=[X1(t),X2(t),X3(t),U(t)]. The inputs to the dynamics file are X and
% U, so we willl take derivatives wrt Y(t).
% First create Derivative inputs for X and U
gX = adigatorCreateDerivInput([Inf,3],struct('vodname','Y','vodsize',[Inf,4],...
  'nzlocs',[(1:3).' (1:3).']));
gU = adigatorCreateDerivInput([Inf,1],struct('vodname','Y','vodsize',[Inf,4],...
  'nzlocs',[1 4]));
% Create Dynamics Vectorized Derivative File
dyn_out = adigator('dynamics',{gX,gU},'dynamics_Y',adigatorOptions('overwrite',1));
% Can extract the sparsity pattern of JF(Y(t)) from dyn_out
I = dyn_out{1}.deriv.nzlocs(:,1); J = dyn_out{1}.deriv.nzlocs(:,2);
time(1) = toc;

% Note, this is done outside of the loop - since its vectorized we only
% need to create this file once.

for i = 1:4
tic;
% ---------------------- Set Up the Problem ----------------------------- %
if i == 1
[probinfo,upperbound,lowerbound,guess,tau] = setupproblem(numintervals(i));
else
% If on second or third mesh, use previous solution as initial guess
[probinfo,upperbound,lowerbound,guess,tau] = ...
  setupproblem(numintervals(i),X,U,fval,tau);
end

% ------------- Project Vectorized Dynamics Non-zero locs --------------- %
N = length(tau);
[II,JJ] = adigatorProjectVectLocs(N,I,J);
% These are not in the native MATLAB order, but rather in the order of the
% outputs from dynamics_Z - need to re-order and store the index so that we
% can use it from within our constraint wrapper.
dFpat = sparse(II,JJ,1:numel(II),3*N,4*N);
% Note: we placed tf at the end of our decision vector so this sparsity
% pattern is also valid for our decision vector which is [X U tf]
[II,JJ,probinfo.dFind] = find(dFpat);

% ----------------- Create Constraint Derivative File ------------------- %
% Here we will be taking the derivative of the file
% vect_cons(X,F,tf,probinfo), so we need to make derivative inputs for each
% of these

% input x has size N by 3 and unrolled jacobian with non-zero derivatives
% along the 1:3*N diagonal elements.
gx = adigatorCreateDerivInput([N 3],...
  struct('vodname','z','vodsize',[length(guess),1],...
  'nzlocs',[(1:3*N).' (1:3*N).']));
% input f has size N by 3 and unrolled jacobian with non-zero derivatives
% defined by [II JJ]
gf = adigatorCreateDerivInput([N 3],...
  struct('vodname','z','vodsize',[length(guess),1],...
  'nzlocs',[II JJ]));
% tf is last element of z
gtf = adigatorCreateDerivInput([1 1],...
  struct('vodname','z','vodsize',[length(guess),1],...
  'nzlocs',[1 4*N+1]));
% Call adigator
adigator('vect_cons',{gx,gf,gtf,probinfo},'vect_cons_z',adigatorOptions('overwrite',1));
% We create a wrapper for this file in vect_conswrap.m
% --------------------------- Call fmincon ------------------------------ %
[z,fval] = ...
  fmincon(@(x)basic_objwrap(x,probinfo),guess,[],[],[],[],...
  lowerbound,upperbound,@(x)vect_conswrap(x,probinfo),options);

% --------------------------- Extract Solution -------------------------- %
t = tau(:).*fval;
X = z(probinfo.xind);
U = z(probinfo.uind);
time(i) = time(i) + toc;

% ---------------------------- Plot Solution ---------------------------- %
figure(i);
subplot(2,1,1)
plot(t,X,'-o');
xlabel('time')
ylabel('states')
legend('x-position','y-position','speed','location','northwest')
title(sprintf(['Mesh %1.0f Solution in ',num2str(time(i)),'s'],i))
subplot(2,1,2);
plot(t,U,'-o');
xlabel('time')
ylabel('control')
end


fprintf(['Total Time Supplying First Derivatives (vectorized): ',...
  num2str(sum(time)),'\n']);