function xdot = TwoLinkSys(t,x)
% Two Link Robot Manipulator System
% Inputs: 
%   x(1) - link 1 position
%   x(2) - link 2 position
%   x(3) - link 1 velocity
%   x(4) - link 2 velocity
%   x(5) - internal filter variable p1
%   x(6) - internal filter variable p2
%   x(7) - stuff we have to numerically integrate for Thetahat1 - z1
%   x(8) - stuff we have to numerically integrate for Thetahat2 - z2
%   x(9) - stuff we have to numerically integrate for Thetahat3 - z3
%   x(10) - stuff we have to numerically integrate for Thetahat4 - z4
%   x(11) - stuff we have to numerically integrate for Thetahat5 - z5
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global probinfo
q1    = x(1);
q2    = x(2);
q1dot = x(3);
q2dot = x(4);
q     = [q1; q2];
qdot  = [q1dot; q2dot];
p     = x(5:6);
z     = x(7:11);

% ROBOT Model from Burg 97
p1  = probinfo.p1;
p2  = probinfo.p2;
p3  = probinfo.p3;
fd1 = probinfo.fd1;
fd2 = probinfo.fd2;

M=zeros(2,2);
M(1,1)=p1+2*p3*cos(q2);
M(1,2)=p2+p3*cos(q2);
M(2,1)=p2+p3*cos(q2);
M(2,2)=p2;

V=zeros(2,2);
V(1,1)=-p3*sin(q2)*q2dot;
V(1,2)=-p3*sin(q2)*(q1dot+q2dot);
V(2,1)= p3*sin(q2)*q1dot;
V(2,2)=0;

Fd= diag([fd1 fd2]);

% Get control law
[u,pdot,zdot] = getDCALcontrol(t,q,p,z);

qdotdot = M\(u-V*qdot-Fd*qdot);

% Define output
xdot = [qdot;qdotdot;pdot;zdot];
end

