% X = SOLVE_CUBIC(a,b,c)
%
% Solve the cubic equation
%
%   x^3 + a*x^2 + b*x + c = 0
%
% with real-valued coefficients a, b, c analytically using Cardano's method.
%
% ARGUMENTS
% a      [n,1]  coefficient vector
%     or [1,1]  coefficient scalar
% b      [n,1]  coefficient vector
%     or [1,1]  coefficient scalar
% c      [n,1]  coefficient vector
%     or [1,1]  coefficient scalar
% 
% RESULTS
% X      [n,3]  result matrix, each column contains a solution; the real one
%               comes first
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 11

function X = solve_cubic(a,b,c)

n = max([numel(a),numel(b),numel(c)]); X = zeros(n,3);         % allocate memory
p = (3*b(:)-a(:).*a(:))/3;             % y^3 + p*y + q = 0, substitute x = y-a/3
q = 2*a(:).*a(:).*a(:)/27 - a(:).*b(:)/3 + c(:);
D = p.*p.*p/27 + q.*q/4; d = sqrt(D);                             % discriminant

id = D(:)>=0;                                         % nonnegative discriminant
u = -q(id)/2+d(id); u = sign(u).*abs(u).^(1/3);
v = -q(id)/2-d(id); v = sign(v).*abs(v).^(1/3);
r1 = (-1+1i*sqrt(3))/2; r2 = (-1-1i*sqrt(3))/2;
X(id,1) = u+v; X(id,2) = r1*u+r2*v;  X(id,3) = r2*u+r1*v;

id = D(:)<0;                                             % negative discriminant
r = sqrt(-p(id).*p(id).*p(id)/27);
phi = acos(-q(id)./(2*r)); al = 2*r.^(1/3);
X(id,1) = al.*cos(phi/3);
X(id,2) = al.*cos(phi/3+2*pi/3);
X(id,3) = al.*cos(phi/3+4*pi/3);

for i=1:3, X(:,i) = X(:,i) - a(:)/3; end                        % backsubstitute