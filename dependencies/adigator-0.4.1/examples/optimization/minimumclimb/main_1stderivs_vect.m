% This function solves the minimum time to climb problem by supplying first
% derivatives which take advantage of the vectorized nature of the problem.
% A derivative file for the dynamics is created as dynamics_yvect and these
% derivatives are used to build the constraint jacobian within the file
% conswrap.m
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
options = optimset('Algorithm','interior-point','MaxFunEvals',50000,...
  'GradObj','on','GradConstr','on','Display','iter','MaxIter',1000);
time = zeros(4,1);
numintervals = [10,20,40];
[probinfo,upperbound,lowerbound,guess,tau] = setupproblem(numintervals(1));
% ------------------- Create Dynamics Derivative File ------------------- %
% Since we are doing this vectorized, we only need to take derivatives once
% Inputs to dynamics file are X,U,probinfo
% X is size N by n
% U is size N by m
% Let y(t) = [x(t) u(t)], Y = [X U], taking derivatives of dynamics file
% wrt y(t).
% Note: our decision vector, z, is of form z = [Y(:);tf], and dynamics file
% is not a function of tf.
tic
n = 3; m = 1; 
gX = adigatorCreateDerivInput([Inf n],...
  struct('vodname','Y','vodsize',[Inf (m+n)],...
  'nzlocs',[(1:n).' (1:n).']));
gU = adigatorCreateDerivInput([Inf m],...
  struct('vodname','Y','vodsize',[Inf (m+n)],...
  'nzlocs',[(1:m).' (n+1:m+n).']));
gOut = adigator('dynamics',{gX,gU,probinfo},'dynamics_yvect',...
  adigatorOptions('overwrite',1));
% We can now extract sparsity pattern of dF(t)/dy(t)
iFy = gOut{1}.deriv.nzlocs(:,1);
jFy = gOut{1}.deriv.nzlocs(:,2);
time(1) = toc;

for i = 1:3
tic;
% ---------------------- Set Up the Problem ----------------------------- %
if i > 1
% If on second or third mesh, use previous solution as initial guess
[probinfo,upperbound,lowerbound,guess,tau] = ...
  setupproblem(numintervals(i),X,U,fval,tau);
end
% ----------------------- Project Deriv Indices ------------------------- %
% We have already created the derivative file, but now that we have fixed
% the vectorized dimension, we can find the unrolled locations of the
% unrolled derivative. We do this here rather than on each derivative
% evaluation to save time.
[probinfo.iFy, probinfo.jFy] = adigatorProjectVectLocs(probinfo.N,iFy,jFy);
probinfo.vectders = 1;

% --------------------------- Call fmincon ------------------------------ %
[z,fval] = ...
  fmincon(@(x)basic_objwrap(x,probinfo),guess,[],[],[],[],...
  lowerbound,upperbound,@(x)conswrap(x,probinfo),options);

% --------------------------- Extract Solution -------------------------- %
t = tau(:).*z(end);
X = z(probinfo.xind);
U = z(probinfo.uind);
time(i) = time(i) + toc;
disp(z(end))
% ---------------------------- Plot Solution ---------------------------- %
figure(i);
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


fprintf(['Total Time Supplying First Derivatives (vectorized): ',...
  num2str(sum(time)),'\n']);