
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
% written in a vectorized manner which computes F(Y(t)), where
% Y(t)=[X1(t),X2(t),X3(t),U(t)]. The inputs to the dynamics file are X and
% U, so we willl take derivatives wrt Y(t).
% First create Derivative inputs for X and U
gX = adigatorCreateDerivInput([Inf,3],struct('vodname','Y','vodsize',[Inf,4],...
  'nzlocs',[(1:3).' (1:3).']));
gU = adigatorCreateDerivInput([Inf,1],struct('vodname','Y','vodsize',[Inf,4],...
  'nzlocs',[1 4]));
% Create Vectorized Dynamics 1st Derivative File
dyn_out1 = adigator('dynamics',{gX,gU},'dynamics_Y',adigatorOptions('overwrite',1));
% Can extract the sparsity pattern of dF(t)/dY(t) from dyn_out
I1 = dyn_out1{1}.deriv.nzlocs(:,1); J1 = dyn_out1{1}.deriv.nzlocs(:,2);
time(1) = toc;
% Create Vectorized Dynamics 2nd Derivative File
gXs.f = gX; gUs.f = gU;
% Need to make aux inputs for derivative fields since we dont know their
% vectorized dimension
gXs.dY = adigatorCreateAuxInput([Inf 3]);
gUs.dY = adigatorCreateAuxInput([Inf 1]);
dyn_out2 = adigator('dynamics_Y',{gXs,gUs},'dynamics_YY',...
  adigatorOptions('overwrite',1,'comments',0));

% Now, we will need the sparsity pattern of F.dY wrt Y(t)
% NOTE: this is a different sparsity pattern than than d^2F/dY^2, it is the
% sparsity pattern of a vector of the non-zero locations of dF/dY wrt Y.
I2 = dyn_out2{1}.dY.deriv.nzlocs(:,1); J2 = dyn_out2{1}.dY.deriv.nzlocs(:,2);
% Note, this is done outside of the loop - since its vectorized we only
% need to create this file once.



for i = 1:4
start = tic;
% ---------------------- Set Up the Problem ----------------------------- %
if i == 1
[probinfo,upperbound,lowerbound,guess,tau] = setupproblem(numintervals(i));
else
% If on second or third mesh, use previous solution as initial guess
[probinfo,upperbound,lowerbound,guess,tau] = ...
  setupproblem(numintervals(i),X,U,fval,tau);
end

% ------ Project Vectorized Dynamics Non-zero locs - 1st Deriv ---------- %
N = length(tau);
[II1,JJ1] = adigatorProjectVectLocs(N,I1,J1);
% These are not in the native MATLAB order, but rather in the order of the
% outputs from dynamics_Z - need to re-order and store the index so that we
% can use it from within our constraint wrapper.
dFpat = sparse(II1,JJ1,1:numel(II1),3*N,4*N);
% Note: we placed tf at the end of our decision vector so this sparsity
% pattern is also valid for our decision vector which is [X U tf]
[II1,JJ1,probinfo.dFind] = find(dFpat);

% ------ Project Vectorized Dynamics Non-zero locs - 2nd Deriv ---------- %
N = length(tau);
[II2,JJ2] = adigatorProjectVectLocs(N,I2,J2);
% These are not in the native MATLAB order, but rather in the order of the
% outputs from dynamics_Y - need to re-order and store the index so that we
% can use it from within our constraint wrapper.
dF2pat = sparse(II2,JJ2,1:numel(II2),numel(II1),4*N);
% Need to re-order this because we re0order the first derivative.
dF2pat = dF2pat(probinfo.dFind,:);
% Note: this is the sparsity pattern of the derivative of the non-zero
% elemnts of dF/dY (i.e. the non-zero elements of dFpat) taken wrt Y - thus
% this unrolled Jacobian is of size numel(II1) by 4*N.
[II2,JJ2,probinfo.dF2ind] = find(dF2pat);

% ---------------- Create 1st Derivative Constraint File ---------------- %
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
  'nzlocs',[II1 JJ1]));
% tf is last element of z
gtf = adigatorCreateDerivInput([1 1],...
  struct('vodname','z','vodsize',[length(guess),1],...
  'nzlocs',[1 4*N+1]));
% Call adigator
Coutput = adigator('vect_cons',{gx,gf,gtf,probinfo},'vect_cons_z',...
  adigatorOptions('overwrite',1));
% We create a wrapper for this file in vect_conswrap.m

% ---------- Take Derivative of Lagrangian Gradient File ---------------- %
% We have built a file which creates the Lagrangian Gradient,
% vect_laggrad(Xf,dX,Ff,dF,tff,dtf,lambda,probinfo), by calling the
% function vect_cons_z. We will take the derivative of this file in order
% to build our Lagrangian Hessian. First, we must create inputs for the
% dX,dF,dtf.
gdx = ones(3*N,1); % dx has no derivative wrt z
gdf = adigatorCreateDerivInput([numel(II1), 1],...
  struct('vodname','z','vodsize',[length(guess),1],...
  'nzlocs',[II2,JJ2])); % df has derivatives wrt z defined by [II2 JJ2]
gdtf = 1; % dtf has no derivative wrt z
glambda = adigatorCreateAuxInput(Coutput{1}.size);
% Call adigator
gH2 = adigator('vect_laggrad',{gx,gdx,gf,gdf,gtf,gdtf,glambda,probinfo},...
  'vect_laggrad_z',adigatorOptions('overwrite',1));

% ------------------------ Define Hessian File -------------------------- %
% We wrote a wrapper for vect_laggrad_z and called it vect_laggradwrap.m -
% this will create the Lagrangian Hessian from vect_laggrad_z.
options = optimset(options,'Hessian','user-supplied','HessFcn',...
  @(z,lambda)vect_laggradwrap(z,lambda,probinfo));

% --------------------------- Call fmincon ------------------------------ %
[z,fval,exitflag,output,lambda] = ...
  fmincon(@(x)basic_objwrap(x,probinfo),guess,[],[],[],[],...
  lowerbound,upperbound,@(x)vect_conswrap(x,probinfo),options);

% --------------------------- Extract Solution -------------------------- %
t = tau(:).*fval;
X = z(probinfo.xind);
U = z(probinfo.uind);
time(i) = time(i) + toc(start);

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

fprintf(['Total Time Supplying Second Derivatives (vectorized): ',...
  num2str(sum(time)),'\n']);