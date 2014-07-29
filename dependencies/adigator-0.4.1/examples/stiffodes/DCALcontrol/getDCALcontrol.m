function [u,pdot,zdot] = getDCALcontrol(t,q,p,z)
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global probinfo
K     = probinfo.DCAL.K;
Gamma = probinfo.DCAL.Gamma;

% Add noise to signal if desired
if probinfo.noiseflag
  noise = interp1(probinfo.time.',probinfo.noise.',t).';
  q = q.*(1+noise);
end

% Get Desired Trajectory plus time derivatives
qd = getqd_dtdtdt(struct('f',t,'dt',1));

% Compute Error 
e  = qd.f-q;

% Build pdot
pdot = -(K+1)*p + (K.^2 + 1)*e;

% Build filtered error, ef
ef   = -K*e+p;

% Q    = [q1 q2 q_1 q_2 q__1 q__2]'
% dQdt = [q_1 q_2 q__1 q__2 q___1 q___2]'
Q.f    = [qd.f; qd.dt; qd.dtdt];
Q.dt   = [qd.dt; qd.dtdt; qd.dtdtdt];

% Build Yd
Yd = getYd_dt(Q);
Yddot = zeros(Yd.dt_size);
Yddot(sub2ind(Yd.dt_size,Yd.dt_location(:,1),Yd.dt_location(:,2))) = Yd.dt;
Yd = Yd.f;

% Build zdot
zdot = Yd.'*(e + ef) - Yddot.'*e;

% Build Thetahat
Thetahat = Gamma*(z + Yd.'*e);

% Define Control Law
u = Yd*Thetahat + -K*ef + e;
end
