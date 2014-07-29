% This function solves the brachistochrone problem by supplying first
% derivatives without taking advantage of the vectorized nature of the
% problem. A wrapper file for the constraint and objective functions (which
% call derivative files) are created and supplied in basic_conswrap.m and
% basic_objwrap.m, respectively.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%options = optimset('Algorithm','interior-point','MaxFunEvals',50000,...
%  'GradObj','on','GradConstr','on','Display','iter');
time = zeros(4,1);
numintervals = [5,10,20,40];
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

% ----------------- Create Constraint Derivative File ------------------- %
% NOTE: fmincon is given the function basic_conswrap which builds the
% Jacobian by calling basic_cons_z
gz = adigatorCreateDerivInput([length(guess), 1],'z');
adigator('basic_cons',{gz,probinfo},'basic_cons_z',adigatorOptions('overwrite',1));
% Also note, we must create a new derivative file each time we change
% meshes

% --------------------------- Call fmincon ------------------------------ %
[z,fval] = ...
  fmincon(@(x)basic_objwrap(x,probinfo),guess,[],[],[],[],...
  lowerbound,upperbound,@(x)basic_conswrap(x,probinfo),options);

% --------------------------- Extract Solution -------------------------- %
t = tau(:).*fval;
X = z(probinfo.xind);
U = z(probinfo.uind);
time(i) = toc;

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


fprintf(['Total Time Supplying First Derivatives (non-vectorized): ',...
  num2str(sum(time)),'\n']);