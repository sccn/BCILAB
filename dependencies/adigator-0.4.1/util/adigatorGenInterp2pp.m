function pp = adigatorGenInterp2pp(X,Y,Z,method)
% function pp = adigatorGenInterp2pp(X,Y,Z,method)
% NON-OVERLOADED VERSION
% This function generates a 2-D piecewise polynomial, given data points
% X,Y, and Z. Inputs to this should be the same as to the MATLAB interp2
% function, with the exception that XI and YI are not to be specified.
% Furthermore, the output of this is similar to the command pp =
% interp1(x,y,method,'pp'), except that the polynomial must be evaluated
% using the command ZI = adigatorEvalInterp2pp(pp,XI,YI).
%
% Methods: 
%          'linear'  - bilinear interpolation
%          'spline'  - spline interpolation
%          'cubic'   - bicubic interpolation as long as the data is
%                      uniformly spaced, otherwise the same as 'spline'
% NOTE: This is not written for the 'nearest' case - the coefficients are
% simply Z.
% 
%
% The following two lines of code:
%
%     pp = adigatorGenInterp2pp(X,Y,Z,method);
%     ZI = adigatorEvalInterp2pp(pp,XI,YI);
%
% Are (roughly) equivalent to the single line of code:
%
%     ZI = interp2(X,Y,Z,XI,YI,method);
%
% Where the difference is due to rounding errors since when using the 
% adigatorGenInterp2pp and adigatorEvalInterp2pp commands, the 2-D
% coefficients are calculated, whereas MATLAB interp2 essentially performs
% interp1 twice, thus generating 1-D coefficients twice.
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
% 
% See also interp2, interp1, ppval, adigatorEvalInterp2pp

[msg, X,Y,Z,~] = xyzchk(X,Y,Z);
error(msg)
pp.form = 'adigatorpp2';
pp.xbreaks = X(1,:);
pp.ybreaks = Y(:,1).';
switch method(1)
  case 'l'
    % Bi-Linear
    pp.coefs  = linearcoefs2(X,Y,Z);
    pp.xorder = 2;
    pp.yorder = 2;
  case 's'
    % Spline
    pp.coefs  = splinecoefs2(X,Y,Z);
    pp.xorder = 4;
    pp.yorder = 4;
  case 'c'
    % Cubic
    xdif = diff(pp.xbreaks); ydif = diff(pp.ybreaks);
    if all(xdif(1)==xdif) && all(ydif(1)==ydif)
      % Data Equally spaced - use cubic
      pp.coefs = cubiccoefs2(X,Y,Z);
    else
      % Data not equally spaced - use spline
      pp.coefs = splinecoefs2(X,Y,Z);
    end
    pp.xorder = 4;
    pp.yorder = 4;
  otherwise
    error('invalid method')
end

end

function D = linearcoefs2(x,y,z)
% Generate bi-linear coefficients

N = size(x,2);
M = size(y,1);

x = x(1:M-1,:);
y = y(:,1:N-1);

difx = diff(x,1,2);
dify = diff(y,1,1);

difz = diff(z,1,2);

a1 = difz(1:M-1,:)./difx;
a2 = z(1:M-1,1:N-1);

b1 = difz(2:M,:)./difx;
b2 = z(2:M,1:N-1);


d1 = (b1-a1)./dify;
d2 = (b2-a2)./dify;
d3 = a1;
d4 = a2;

%D = zeros(M-1,N-1,2,2);
% D(:,:,1,1) = d1;
% D(:,:,1,2) = d2;
% D(:,:,2,1) = d3;
% D(:,:,2,2) = d4;
% z = D(1,1)*X^1*Y^1 + D(1,2)*X^0*Y^1
%   + D(2,1)*X^1*Y^0 + D(2,2)*X^0*Y^0

D = cell(2,2);
D{1,1} = d1;
D{1,2} = d2;
D{2,1} = d3;
D{2,2} = d4;

end

function coefs = cubiccoefs2(x,y,z)
%   Based on "Cubic Convolution Interpolation for Digital Image
%   Processing", Robert G. Keys, IEEE Trans. on Acoustics, Speech, and
%   Signal Processing, Vol. 29, No. 6, Dec. 1981, pp. 1153-1160.

[M,N] = size(z);

% x and y must have equally spaced intervals, or else MATLAB calls spline2
dx = (x(1,N)-x(1,1))/(N-1);
dy = (y(M,1)-y(1,1))/(M-1);

% create boundary conditions for Z matrix
Zc = [3*z(:,1)-3*z(:,2)+z(:,3), z, 3*z(:,N)-3*z(:,N-1)+z(:,N-2)];

% Zc is M x N+2
% define indices, such that evaluating on interval [xi xj]
H = 1:N-1;
I = 2:N;
J = 3:N+1;
K = 4:N+2;
% create A, B, C, D matrices, R = A(x-xi)^3+B(x-xi)^2+C(x-xi)+D
A = ( -Zc(:,H) + 3*Zc(:,I) - 3*Zc(:,J) + Zc(:,K) )/(2*dx^3);
B = ( 2*Zc(:,H) - 5*Zc(:,I) + 4*Zc(:,J) - Zc(:,K) )/(2*dx^2);
C = ( -Zc(:,H) + Zc(:,J) )/(2*dx);
D = Zc(:,I);
% A,B,C,D are M x N-1

% create boundary conditions for R coefficients
Ac = [3*A(1,:)- 3*A(2,:)+A(3,:); A; 3*A(M,:)-3*A(M-1,:)+A(M-2,:)];
Bc = [3*B(1,:)- 3*B(2,:)+B(3,:); B; 3*B(M,:)-3*B(M-1,:)+B(M-2,:)];
Cc = [3*C(1,:)- 3*C(2,:)+C(3,:); C; 3*C(M,:)-3*C(M-1,:)+C(M-2,:)];
Dc = [3*D(1,:)- 3*D(2,:)+D(3,:); D; 3*D(M,:)-3*D(M-1,:)+D(M-2,:)];
% Ac,Bc,Cc,Dc are M+2 x N-1
% redefine indices, such that evaluating on interval [yi yj]
H = 1:M-1;
I = 2:M;
J = 3:M+1;
K = 4:M+2;

% Define output coefficients Z = Ay*(y-yi)^3 + By(y-yi)^2 + Cy(y-yi) + Di, 
% where Ay, By, Cy, Dy are cubic in (x-xi);
% coefs = zeros(M-1,N-1,4,4);
% % from Ay
% coefs(:,:,1,1) = ( -Ac(H,:) + 3*Ac(I,:) - 3*Ac(J,:) + Ac(K,:) )/(2*dy^3);
% coefs(:,:,1,2) = ( -Bc(H,:) + 3*Bc(I,:) - 3*Bc(J,:) + Bc(K,:) )/(2*dy^3);
% coefs(:,:,1,3) = ( -Cc(H,:) + 3*Cc(I,:) - 3*Cc(J,:) + Cc(K,:) )/(2*dy^3);
% coefs(:,:,1,4) = ( -Dc(H,:) + 3*Dc(I,:) - 3*Dc(J,:) + Dc(K,:) )/(2*dy^3);
% % from By
% coefs(:,:,2,1) = ( 2*Ac(H,:) - 5*Ac(I,:) + 4*Ac(J,:) - Ac(K,:) )/(2*dy^2);
% coefs(:,:,2,2) = ( 2*Bc(H,:) - 5*Bc(I,:) + 4*Bc(J,:) - Bc(K,:) )/(2*dy^2);
% coefs(:,:,2,3) = ( 2*Cc(H,:) - 5*Cc(I,:) + 4*Cc(J,:) - Cc(K,:) )/(2*dy^2);
% coefs(:,:,2,4) = ( 2*Dc(H,:) - 5*Dc(I,:) + 4*Dc(J,:) - Dc(K,:) )/(2*dy^2);
% % from Cy
% coefs(:,:,3,1) = ( -Ac(H,:) + Ac(J,:) )/(2*dy);
% coefs(:,:,3,2) = ( -Bc(H,:) + Bc(J,:) )/(2*dy);
% coefs(:,:,3,3) = ( -Cc(H,:) + Cc(J,:) )/(2*dy);
% coefs(:,:,3,4) = ( -Dc(H,:) + Dc(J,:) )/(2*dy);
% % from Dy
% coefs(:,:,4,1) = Ac(I,:);
% coefs(:,:,4,2) = Bc(I,:);
% coefs(:,:,4,3) = Cc(I,:);
% coefs(:,:,4,4) = Dc(I,:);

coefs = cell(4,4);
% from Ay
coefs{1,1} = ( -Ac(H,:) + 3*Ac(I,:) - 3*Ac(J,:) + Ac(K,:) )/(2*dy^3);
coefs{1,2} = ( -Bc(H,:) + 3*Bc(I,:) - 3*Bc(J,:) + Bc(K,:) )/(2*dy^3);
coefs{1,3} = ( -Cc(H,:) + 3*Cc(I,:) - 3*Cc(J,:) + Cc(K,:) )/(2*dy^3);
coefs{1,4} = ( -Dc(H,:) + 3*Dc(I,:) - 3*Dc(J,:) + Dc(K,:) )/(2*dy^3);
% from By
coefs{2,1} = ( 2*Ac(H,:) - 5*Ac(I,:) + 4*Ac(J,:) - Ac(K,:) )/(2*dy^2);
coefs{2,2} = ( 2*Bc(H,:) - 5*Bc(I,:) + 4*Bc(J,:) - Bc(K,:) )/(2*dy^2);
coefs{2,3} = ( 2*Cc(H,:) - 5*Cc(I,:) + 4*Cc(J,:) - Cc(K,:) )/(2*dy^2);
coefs{2,4} = ( 2*Dc(H,:) - 5*Dc(I,:) + 4*Dc(J,:) - Dc(K,:) )/(2*dy^2);
% from Cy
coefs{3,1} = ( -Ac(H,:) + Ac(J,:) )/(2*dy);
coefs{3,2} = ( -Bc(H,:) + Bc(J,:) )/(2*dy);
coefs{3,3} = ( -Cc(H,:) + Cc(J,:) )/(2*dy);
coefs{3,4} = ( -Dc(H,:) + Dc(J,:) )/(2*dy);
% from Dy
coefs{4,1} = Ac(I,:);
coefs{4,2} = Bc(I,:);
coefs{4,3} = Cc(I,:);
coefs{4,4} = Dc(I,:);
end

function coefs = splinecoefs2(x,y,z)

[M,N] = size(z);

dx    = diff(x,1,2);
dy    = diff(y,1,1);
difzx = diff(z,1,2);
difzy = diff(z,1,1);

dzdx   = spline2slopes(x(1,:),z,2);
dzdy   = spline2slopes(y(:,1),z,1);
dzdxdy = spline2slopes(y(:,1),dzdx,1);

dx = dx(1:M-1,:);
dy = dy(:,1:N-1);

A1 = difzx(1:M-1,:)./dx;
B1 = (A1 - dzdx(1:M-1,1:N-1))./dx;
C1 = (dzdx(1:M-1,2:N) - A1)./dx;
D1 = (C1 - B1)./dx;
E1 = (2.*B1 - C1);

A2 = difzx(2:M,:)./dx;
B2 = (A2 - dzdx(2:M,1:N-1))./dx;
C2 = (dzdx(2:M,2:N) - A2)./dx;
D2 = (C2 - B2)./dx;
E2 = (2.*B2 - C2);

A1y = (dzdy(1:M-1,2:N)- dzdy(1:M-1,1:N-1))./dx;
B1y = (A1y - dzdxdy(1:M-1,1:N-1))./dx;
C1y = (dzdxdy(1:M-1,2:N) - A1y)./dx;
D1y = (C1y - B1y)./dx;
E1y = (2.*B1y - C1y);

A2y = (dzdy(2:M,2:N)-dzdy(2:M,1:N-1))./dx;
B2y = (A2y - dzdxdy(2:M,1:N-1))./dx;
C2y = (dzdxdy(2:M,2:N) - A2y)./dx;
D2y = (C2y - B2y)./dx;
E2y = (2.*B2y - C2y);

% Build Coefs--------------------------------------------------------------
% coefs = zeros(M-1,N-1,4,4);
% % from D3
% coefs(:,:,1,1) = (-2*(D2 - D1)./dy + D2y + D1y)./(dy.^2);
% coefs(:,:,2,1) = (-2*(E2-E1)./dy + E2y + E1y)./(dy.^2);
% coefs(:,:,3,1) = (-2*(dzdx(2:M,1:N-1) - dzdx(1:M-1,1:N-1))./dy +...
%     dzdxdy(2:M,1:N-1) + dzdxdy(1:M-1,1:N-1))./(dy.^2);
% coefs(:,:,4,1) = (-2*difzy(:,1:N-1)./dy + dzdy(2:M,1:N-1) +...
%     dzdy(1:M-1,1:N-1))./(dy.^2);
% % from E3
% coefs(:,:,1,2) = (3*(D2-D1)./dy - 2*D1y - D2y)./dy;
% coefs(:,:,2,2) = (3*(E2-E1)./dy - 2*E1y - E2y)./dy;
% coefs(:,:,3,2) = (3*(dzdx(2:M,1:N-1) - dzdx(1:M-1,1:N-1))./dy -...
%     2*dzdxdy(1:M-1,1:N-1) - dzdxdy(2:M,1:N-1))./dy;
% coefs(:,:,4,2) = (3*difzy(:,1:N-1)./dy - 2*dzdy(1:M-1,1:N-1) - dzdy(2:M,1:N-1))./dy;
% % from R1y
% coefs(:,:,1,3)  = D1y;
% coefs(:,:,2,3) = E1y;
% coefs(:,:,3,3) = dzdxdy(1:M-1,1:N-1);
% coefs(:,:,4,3) = dzdy(1:M-1,1:N-1);
% % from R1
% coefs(:,:,1,4) = D1;
% coefs(:,:,2,4) = E1;
% coefs(:,:,3,4) = dzdx(1:M-1,1:N-1);
% coefs(:,:,4,4) = z(1:M-1,1:N-1);
% 
% % from D3
% coefs(:,:,1,1) = (-2*(D2 - D1)./dy + D2y + D1y)./(dy.^2);
% coefs(:,:,1,2) = (-2*(E2-E1)./dy + E2y + E1y)./(dy.^2);
% coefs(:,:,1,3) = (-2*(dzdx(2:M,1:N-1) - dzdx(1:M-1,1:N-1))./dy +...
%     dzdxdy(2:M,1:N-1) + dzdxdy(1:M-1,1:N-1))./(dy.^2);
% coefs(:,:,1,4) = (-2*difzy(:,1:N-1)./dy + dzdy(2:M,1:N-1) +...
%     dzdy(1:M-1,1:N-1))./(dy.^2);
% % from E3
% coefs(:,:,2,1) = (3*(D2-D1)./dy - 2*D1y - D2y)./dy;
% coefs(:,:,2,2) = (3*(E2-E1)./dy - 2*E1y - E2y)./dy;
% coefs(:,:,2,3) = (3*(dzdx(2:M,1:N-1) - dzdx(1:M-1,1:N-1))./dy -...
%     2*dzdxdy(1:M-1,1:N-1) - dzdxdy(2:M,1:N-1))./dy;
% coefs(:,:,2,4) = (3*difzy(:,1:N-1)./dy - 2*dzdy(1:M-1,1:N-1) - dzdy(2:M,1:N-1))./dy;
% % from R1y
% coefs(:,:,3,1)  = D1y;
% coefs(:,:,3,2) = E1y;
% coefs(:,:,3,3) = dzdxdy(1:M-1,1:N-1);
% coefs(:,:,3,4) = dzdy(1:M-1,1:N-1);
% % from R1
% coefs(:,:,4,1) = D1;
% coefs(:,:,4,2) = E1;
% coefs(:,:,4,3) = dzdx(1:M-1,1:N-1);
% coefs(:,:,4,4) = z(1:M-1,1:N-1);

coefs = cell(4,4);
% from D3
coefs{1,1} = (-2*(D2 - D1)./dy + D2y + D1y)./(dy.^2);
coefs{2,1} = (-2*(E2-E1)./dy + E2y + E1y)./(dy.^2);
coefs{3,1} = (-2*(dzdx(2:M,1:N-1) - dzdx(1:M-1,1:N-1))./dy +...
    dzdxdy(2:M,1:N-1) + dzdxdy(1:M-1,1:N-1))./(dy.^2);
coefs{4,1} = (-2*difzy(:,1:N-1)./dy + dzdy(2:M,1:N-1) +...
    dzdy(1:M-1,1:N-1))./(dy.^2);
% from E3
coefs{1,2} = (3*(D2-D1)./dy - 2*D1y - D2y)./dy;
coefs{2,2} = (3*(E2-E1)./dy - 2*E1y - E2y)./dy;
coefs{3,2} = (3*(dzdx(2:M,1:N-1) - dzdx(1:M-1,1:N-1))./dy -...
    2*dzdxdy(1:M-1,1:N-1) - dzdxdy(2:M,1:N-1))./dy;
coefs{4,2} = (3*difzy(:,1:N-1)./dy - 2*dzdy(1:M-1,1:N-1) - dzdy(2:M,1:N-1))./dy;
% from R1y
coefs{1,3}  = D1y;
coefs{2,3} = E1y;
coefs{3,3} = dzdxdy(1:M-1,1:N-1);
coefs{4,3} = dzdy(1:M-1,1:N-1);
% from R1
coefs{1,4} = D1;
coefs{2,4} = E1;
coefs{3,4} = dzdx(1:M-1,1:N-1);
coefs{4,4} = z(1:M-1,1:N-1);

% from D3
coefs{1,1} = (-2*(D2 - D1)./dy + D2y + D1y)./(dy.^2);
coefs{1,2} = (-2*(E2-E1)./dy + E2y + E1y)./(dy.^2);
coefs{1,3} = (-2*(dzdx(2:M,1:N-1) - dzdx(1:M-1,1:N-1))./dy +...
    dzdxdy(2:M,1:N-1) + dzdxdy(1:M-1,1:N-1))./(dy.^2);
coefs{1,4} = (-2*difzy(:,1:N-1)./dy + dzdy(2:M,1:N-1) +...
    dzdy(1:M-1,1:N-1))./(dy.^2);
% from E3
coefs{2,1} = (3*(D2-D1)./dy - 2*D1y - D2y)./dy;
coefs{2,2} = (3*(E2-E1)./dy - 2*E1y - E2y)./dy;
coefs{2,3} = (3*(dzdx(2:M,1:N-1) - dzdx(1:M-1,1:N-1))./dy -...
    2*dzdxdy(1:M-1,1:N-1) - dzdxdy(2:M,1:N-1))./dy;
coefs{2,4} = (3*difzy(:,1:N-1)./dy - 2*dzdy(1:M-1,1:N-1) - dzdy(2:M,1:N-1))./dy;
% from R1y
coefs{3,1}  = D1y;
coefs{3,2} = E1y;
coefs{3,3} = dzdxdy(1:M-1,1:N-1);
coefs{3,4} = dzdy(1:M-1,1:N-1);
% from R1
coefs{4,1} = D1;
coefs{4,2} = E1;
coefs{4,3} = dzdx(1:M-1,1:N-1);
coefs{4,4} = z(1:M-1,1:N-1);
end

function s = spline2slopes(x,z,dim)

if dim == 1
  z = z.';
  x = x.';
end

% Check that data are acceptable and, if not, try to adjust them appropriately
[x,z,sizez,endslopes] = chckxy(x,z);
n = length(x); zd = prod(sizez);

% Generate the cubic spline interpolant in ppform

dd = ones(zd,1); dx = diff(x); divdif = diff(z,[],2)./dx(dd,:);

b=zeros(zd,n);
b(:,2:n-1)=3*(dx(dd,2:n-1).*divdif(:,1:n-2)+dx(dd,1:n-2).*divdif(:,2:n-1));
if isempty(endslopes)
  x31=x(3)-x(1);xn=x(n)-x(n-2);
  b(:,1)=((dx(1)+2*x31)*dx(2)*divdif(:,1)+dx(1)^2*divdif(:,2))/x31;
  b(:,n)=...
    (dx(n-1)^2*divdif(:,n-2)+(2*xn+dx(n-1))*dx(n-2)*divdif(:,n-1))/xn;
else
  x31 = 0; xn = 0; b(:,[1 n]) = dx(dd,[2 n-2]).*endslopes;
end
dxt = dx(:);
c = spdiags([ [x31;dxt(1:n-2);0] ...
  [dxt(2);2*(dxt(2:n-1)+dxt(1:n-2));dxt(n-2)] ...
  [0;dxt(2:n-1);xn] ],[-1 0 1],n,n);

% sparse linear equation solution for the slopes
mmdflag = spparms('autommd');
spparms('autommd',0);
s=b/c;
spparms('autommd',mmdflag);

if dim == 1
  s = s.';
end
end

function [x,y,sizey,endslopes] = chckxy(x,y)
%CHCKXY check and adjust input for SPLINE and PCHIP
%   [X,Y,SIZEY] = CHCKXY(X,Y) checks the data sites X and corresponding data
%   values Y, making certain that there are exactly as many sites as values,
%   that no two data sites are the same, removing any data points that involve 
%   NaNs, reordering the sites if necessary to ensure that X is a strictly
%   increasing row vector and reordering the data values correspondingly,
%   and reshaping Y if necessary to make sure that it is a matrix, with Y(:,j)
%   the data value corresponding to the data site X(j), and with SIZEY the
%   actual dimensions of the given values. 
%   This call to CHCKXY is suitable for PCHIP.
%
%   [X,Y,SIZEY,ENDSLOPES] = CHCKXY(X,Y) also considers the possibility that
%   there are two more data values than there are data sites.
%   If there are, then the first and the last data value are removed from Y
%   and returned separately as ENDSLOPES. Otherwise, an empty ENDSLOPES is
%   returned.  This call to CHCKXY is suitable for SPLINE.
%
%   See also PCHIP, SPLINE.

%   Copyright 1984-2003 The MathWorks, Inc.

% make sure X is a vector:
if length(find(size(x)>1))>1 
  error('MATLAB:chckxy:XNotVector','X must be a vector.') 
end

% ensure X is real
if any(~isreal(x)) 
  error('MATLAB:chckxy:XComplex','The X vector should have real elements.') 
end

% deal with NaN's among the sites:
nanx = find(isnan(x));
if ~isempty(nanx)
   x(nanx) = [];
   warning('MATLAB:chckxy:nan','All data points with NaN as their site will be ignored.')
end

n=length(x);
if n<2 
  error('MATLAB:chckxy:NotEnoughPts','There should be at least two data points.') 
end

% re-sort, if needed, to ensure strictly increasing site sequence:
x=x(:).'; 
dx = diff(x);

if any(dx<0), [x,ind] = sort(x); dx = diff(x); else ind=1:n; end

if ~all(dx), error('MATLAB:chckxy:RepeatedSites','The data sites should be distinct.'), end

% if Y is ND, reshape it to a matrix by combining all dimensions but the last:
sizey = size(y);


while length(sizey)>2&&sizey(end)==1, sizey(end) = []; end


yn = sizey(end); 
sizey(end)=[]; 
yd = prod(sizey);

if length(sizey)>1
   y = reshape(y,yd,yn);
else
   % if Y happens to be a column matrix, change it to the expected row matrix.
   if yn==1
       yn = yd;
       y = reshape(y,1,yn); 
       yd = 1; 
       sizey = yd;
   end
end

% determine whether not-a-knot or clamped end conditions are to be used:
nstart = n+length(nanx);
if yn==nstart
   endslopes = [];
elseif nargout==4&&yn==nstart+2
   endslopes = y(:,[1 n+2]); y(:,[1 n+2])=[];
   if any(isnan(endslopes))
      error('MATLAB:chckxy:EndslopeNaN','The endslopes cannot be NaN.')
   end
   if any(isinf(endslopes))
       error('MATLAB:chckxy:EndslopeInf','The endslopes cannot be Inf.')
   end
else
   error('MATLAB:chckxy:NumSitesMismatchValues',...
        ['The number of sites, ' int2str(nstart), ...
        ', is incompatible with the number of values, ' int2str(yn) '.'])
end

% deal with NaN's among the values:
if ~isempty(nanx)
    y(:,nanx) = [];
end

y=y(:,ind);
nany = find(sum(isnan(y),1));
if ~isempty(nany)
   y(:,nany) = []; x(nany) = [];
   warning('MATLAB:chckxy:IgnoreNaN','All data points with NaN in their value will be ignored.')
   n = length(x);
   if n<2 
     error('MATLAB:chckxy:NotEnoughPts', 'There should be at least two data points.') 
   end
end

end