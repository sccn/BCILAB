function daeout = dynamics(x,u,probinfo)

CONSTANTS = probinfo.CONSTANTS;

CoF  = CONSTANTS.CoF;

h   = x(:,1);
v   = x(:,2);
fpa = x(:,3);

c1  = 392.4;    c2  = 16818;      c3 = 86.138;
c4  = 288.14;   c5  = 6.49;       c6  = 4.0519e9;
c7  = 288.08;   c8  = 5.256;      c9  = 216.64;
c10 = 9.06e8;   c11 = 1.73;       c12 = 0.157;
c13 = 6e-5;     c14 = 4.00936;    c15 = 2.2;


T = zeros(size(h));
p    = zeros(size(h));

ihlow     = h < 11000;
T(ihlow)  = c4 - c5*h(ihlow);
p(ihlow)  = c6*(T(ihlow)./c7).^c8;

ihhigh    = ~ihlow;
T(ihhigh) = c9;
p(ihhigh) = c10* exp(c11 - c12*h(ihhigh));

rho = c3*p./T;
q = 0.5.*rho.*v.*v.*c13;

a = c14.*sqrt(T);       M = v./a;   
M0 = M.^0;  M1 = M.^1;  M2 = M.^2;
M3 = M.^3;  M4 = M.^4;  M5 = M.^5;

numeratorCD0   = CoF(1,1).*M0+CoF(1,2).*M1+CoF(1,3).*M2+CoF(1,4).*M3+CoF(1,5).*M4;
denominatorCD0 = CoF(2,1).*M0+CoF(2,2).*M1+CoF(2,3).*M2+CoF(2,4).*M3+CoF(2,5).*M4;
Cd0            = numeratorCD0./denominatorCD0;

numeratorK   = CoF(3,1).*M0+CoF(3,2).*M1+CoF(3,3).*M2+CoF(3,4).*M3+CoF(3,5).*M4;
denominatorK = CoF(4,1).*M0+CoF(4,2).*M1+CoF(4,3).*M2+CoF(4,4).*M3+CoF(4,5).*M4+CoF(4,6).*M5;
K            = numeratorK./denominatorK;
FD           = q.*(Cd0+K.*((c2.^2).*(c1.^2)./(q.^2)).*(u.^2));

e0 = CoF(5,1).*M0+CoF(6,1).*M1+CoF(7,1).*M2+CoF(8,1).*M3+CoF(9,1).*M4+CoF(10,1).*M5;
e1 = CoF(5,2).*M0+CoF(6,2).*M1+CoF(7,2).*M2+CoF(8,2).*M3+CoF(9,2).*M4+CoF(10,2).*M5;
e2 = CoF(5,3).*M0+CoF(6,3).*M1+CoF(7,3).*M2+CoF(8,3).*M3+CoF(9,3).*M4+CoF(10,3).*M5;
e3 = CoF(5,4).*M0+CoF(6,4).*M1+CoF(7,4).*M2+CoF(8,4).*M3+CoF(9,4).*M4+CoF(10,4).*M5;
e4 = CoF(5,5).*M0+CoF(6,5).*M1+CoF(7,5).*M2+CoF(8,5).*M3+CoF(9,5).*M4+CoF(10,5).*M5;
e5 = CoF(5,6).*M0+CoF(6,6).*M1+CoF(7,6).*M2+CoF(8,6).*M3+CoF(9,6).*M4+CoF(10,6).*M5;
FT = (e0.*h.^0+e1.*h.^1+e2.*h.^2+e3.*h.^3+e4.*h.^4+e5.*h.^5).*c1/c15;

hdot   = v.*sin(fpa);
vdot   = (FT-FD)./c2 - c1.*sin(fpa); 
fpadot = c1.*(u-cos(fpa))./v;   

daeout = [hdot vdot fpadot];