% X = SOLVE_QUARTIC(a,b,c,d)
%
% Solve the quartic equation
%
%   x^4 + a*x^3 + b*x^2 + c*x + d = 0
%
% with real-valued coefficients a, b, c, d analytically using Ferrari's method.
%
% ARGUMENTS
% a      [n,1]  coefficient vector
%     or [1,1]  coefficient scalar
% b      [n,1]  coefficient vector
%     or [1,1]  coefficient scalar
% c      [n,1]  coefficient vector
%     or [1,1]  coefficient scalar
% d      [n,1]  coefficient vector
%     or [1,1]  coefficient scalar
% 
% RESULTS
% X      [n,4]  result matrix, each column contains a solution
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 11

function X = solve_quartic(a,b,c,d)

n = max([numel(a),numel(b),numel(c),numel(d)]); X = zeros(n,4);% allocate memory
p = b(:) - a(:).*a(:)*(3/8);   % y^4 + p*y^2 + r*y + q = 0, substitute x = y-a/4
q = c(:) - a(:).*b(:)/2 + a(:).*a(:).*a(:)/8;
r = d(:) - a(:).*c(:)/4 + a(:).*a(:).*b(:)/16 - a(:).*a(:).*a(:).*a(:)*(3/256);
Z = solve_cubic(2*p, p.*p-4*r, -q.*q);                               % resolvent
Z = sort(Z,2,'descend');                  % bring smallest solution to the right
id = min(abs(Z(:,1:2)),[],2)>1e-12;                            % nonzero entries
W = sqrt(Z); W(id,3) = -q(id)./(W(id,1).*W(id,2));
X = W * [1,1,-1,-1; 1,-1,1,-1; 1,-1,-1,1];
for i=1:4, X(:,i) = X(:,i)/2 - a(:)/4; end                      % backsubstitute