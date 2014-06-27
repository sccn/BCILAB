% X = SOLVE_CUBIC(a,b,c)
%
% Solve the cubic equation
%
%   x^3 + a*x^2 + b*x + c = 0
%
% with real-valued coefficients a, b, c analytically using Cardano's method,
% where first a depressed cubic y^3 - 3*Q*y + 2*R = 0 is obtained by x = y-a/3.
% Notation taken from the numerical recipies in C, .
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
% X      [n,3]  result matrix, each column contains a solution; the first one
%               is real
%
% (c) by Hannes Nickisch, Philips Research, 2012 October 13

function X = solve_cubic(a,b,c)

n = max([numel(a),numel(b),numel(c)]); X = zeros(n,3);         % allocate memory
Q = a(:).*a(:)/9-b(:)/3;
R = a(:).*a(:).*a(:)/27 - a(:).*b(:)/6 + c(:)/2;
R2 = R.*R; Q3 = Q.*Q.*Q; D = R2-Q3; d = sqrt(D);         % D is the discriminant

id = imag(R)==0 & imag(Q)==0 & D<0;                   % we have three real roots
th = acos(R(id)./sqrt(Q3(id))); al = -2*sqrt(Q(id));
X(id,1) = al.*cos( th      /3);
X(id,2) = al.*cos((th+2*pi)/3);
X(id,3) = al.*cos((th-2*pi)/3);

id = ~id;                                      % we do not have three real roots
s = 2*(real(R(id))>0)-1; A = -s.*(s.*R(id)+d(id)).^(1/3);    % have correct sign
B = Q(id)./A; A(abs(A)<1e-10) = 0; B(abs(A)<1e-10) = 0;

X(id,1) =   A+B;
X(id,2) = -(A+B + 1i*sqrt(3)*(A-B))/2;
X(id,3) = -(A+B - 1i*sqrt(3)*(A-B))/2;

for i=1:3, X(:,i) = X(:,i) - a(:)/3; end                        % backsubstitute