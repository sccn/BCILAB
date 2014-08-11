% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
close all
global probinfo

% -------------------- Generate Some Derivative Files ------------------- %
% NOTE: time derivatives of qd and partial derivatives of Yd are generated
% using automatic differentiation software I wrote, the command:
% adigator(function,input size,derivative file name)
% generates a derivative file withe specified name that takes the partial
% derivative of the output wrt the input of the defined function file.

% Need 3rd time derivative of qd
fprintf('\nGenerating third time derivative of qd..\n');
opts = adigatorOptions('overwrite',1,'echo',0);
opts2 = opts;
opts2.comments = 0;
t = adigatorCreateDerivInput([1 1],'t');
adigator('getqd',{t},'getqd_dt',opts);
t = struct('f',t,'dt',1);
adigator('getqd_dt',{t},'getqd_dtdt',opts2);
adigator('getqd_dtdt',{t},'getqd_dtdtdt',opts2);
% Need time derivative of Yd(Q) where Q = [q1 q2 q_1 q_2 q__1 q__2]
fprintf('Generating time derivative of Yd..\n');
Q = adigatorCreateDerivInput([6 1],struct('vodname','t','vodsize',[1 1],...
  'nzlocs',[(1:6).' ones(6,1)]));
adigator('getYd',{Q},'getYd_dt',opts);


% ----------------------- Define Running Parameters --------------------- %
supplyderiv = 1; % Set to 0 to not supply adigator derivs
displayplots = 1;
tf = 50;
h = .2; % 5 points per second
t = 0:.2:tf;
N = length(t);

t = linspace(0,tf,N);
qd = getqd(t);
qd = qd.';
probinfo.noiseflag = 0;       % 0 => no noise, 1 => noise
% Cant put random variables into ode or will never be able to integrate,
% we define a random noise prior to the simulation and interpolate off of
% this to simulate noise.
probinfo.noise = (rand(2,N) - rand(2,N)).*.05;
probinfo.time  = t;

fprintf('\nRunning simulation with following parameters:\n');
fprintf('Noise = %1.0f\n',probinfo.noiseflag);
fprintf('Final time = %1.0f\n',tf);
fprintf('Number of points = %1.0f\n',N);
fprintf('Supplying Derivatives to ODE solver: %1.0f\n',supplyderiv);
fprintf('Display plots: %1.0f\n\n',displayplots);

% -------------------------- Set Robot Params --------------------------- %
probinfo.p1  = 3.473;
probinfo.p2  = 0.193;
probinfo.p3  = 0.242;
probinfo.fd1 = 5.3;
probinfo.fd2 = 1.1;
theta = [probinfo.p1; probinfo.p2; probinfo.p3;...
  probinfo.fd1; probinfo.fd2];

% -------------------------- Set DCAL Params ---------------------------- %
probinfo.DCAL.Gamma  = diag([5 5 30 80 10]);
probinfo.DCAL.K      = 250;


% ------------------------- Set Initial Conditions ---------------------- %
Q0 = zeros(4,1); % System starts from rest
p0 = zeros(2,1);
Thetahat0 = .5*theta;
z0 = probinfo.DCAL.Gamma\Thetahat0;
X0 = [Q0;p0;z0];

% ------------------ Generate Main Derivative File ---------------------- %
if supplyderiv
  fprintf('Generating Derivatives of ODE..\n\n');
  xg = adigatorCreateDerivInput([11 1],'x');
  tg = adigatorCreateAuxInput([1 1]);
  gout = adigator('TwoLinkSys',{tg,xg},'TwoLinkSys_x',opts);
  S = sparse(gout{1}.deriv.nzlocs(:,1),gout{1}.deriv.nzlocs(:,2),...
    ones(size(gout{1}.deriv.nzlocs(:,1),1),1),11,11);
  options = odeset('Jacobian',@adigatorodewrapper,'Jpattern',S,'Stats','on');
else
  options = odeset('Stats','on');
end
% ------------------------- Run Simulation ------------------------------ %
fprintf('Calling ODE solver..\n\n')
tic
[t,X] = ode15s(@TwoLinkSys,t,X0,options);
runtime = toc;
fprintf(['Simulation time = ',num2str(runtime),'\n']);

% --------------------------- Extract Outputs --------------------------- %
q  = X(:,1:2);
q_ = X(:,3:4);
p  = X(:,5:6);
z  = X(:,7:11);

% --------------------- Build Control and Thetahat ---------------------- %
e = (qd-q);
if probinfo.noiseflag
   enoise = qd - q.*(1+probinfo.noise.'); 
end
Gamma = probinfo.DCAL.Gamma;
K     = probinfo.DCAL.K;
if probinfo.noiseflag
  ef    = -K*enoise + p;
else
  ef    = -K*e + p;
end
control  = zeros(N,2);
thetahat = zeros(N,5);
for i = 1:N
  if probinfo.noiseflag
    ei = enoise(i,:).';
  else
    ei = e(i,:).';
  end
  zi = z(i,:).';
  
  % Get Q
  qds = getqd_dtdt(struct('f',t(i),'dt',1));
  Q = [qds.f; qds.dt; qds.dtdt];
  
  % Get Yd
  Yd = getYd(Q);
  
  % Calculate Thetahat
  Thetahati = Gamma*(zi + Yd.'*ei);
  thetahat(i,:) = Thetahati;
  
  % Calculate Control
  efi = ef(i,:).';
  control(i,:) = Yd*Thetahati - K*efi + ei;
end
theta = repmat(theta.',[N 1]);

% --------------------------- Generate Plots ---------------------------- %
if displayplots
figure(1);
plot(t,e(:,1)*180./pi,'-');
xlabel('time (s)');
ylabel('link1 error (deg)');

figure(2)
plot(t,e(:,2)*180./pi,'-');
xlabel('time (s)');
ylabel('link2 error (deg)');

figure(3)
plot(t,control(:,1),'-');
xlabel('time (s)');
ylabel('\tau_1 (Nm)');

figure(4)
plot(t,control(:,2),'-');
xlabel('time (s)');
ylabel('\tau_2 (Nm)');

figure(5)
plot(t,theta(:,1),'b-',t,thetahat(:,1),'b--',...
  t,theta(:,4),'g-',t,thetahat(:,4),'g--');
xlabel('time (s)');
ylabel('adaptive estimates');
legend('p_1','p_1 est','f_d_1','f_d_1 est','Location','Best');

figure(6)
plot(t,theta(:,2),'b-',t,thetahat(:,2),'b--',...
  t,theta(:,3),'g-',t,thetahat(:,3),'g--',t,theta(:,5),'r-',t,thetahat(:,5),'r--');
xlabel('time (s)');
ylabel('adaptive estimates');
legend('p_2','p_2 est','p_3','p_3 est','f_d_2','f_d_2 est','Location','Best');
end