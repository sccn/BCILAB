function [x1,beta,v] = gen_hh(x)
%GEN_HH Generate a Householder transformation.
%
% [x1,beta,v] = gen_hh(x)
%
% Given a vector x, gen_hh computes the scalar beta and the vector v
% determining a Householder transformation
%    H = (I - beta*v*v'),
% such that H*x = +-norm(x)*e_1. x1 is the first element of H*x.

% Based on routine 'house' from Nicholas J. Higham's Test Matrix Toolbox.
% Per Christian Hansen, IMM, Sept. 13, 2001.

alpha = norm(x)*(sign(x(1)) + (x(1)==0));
v = x;
if (alpha==0),
  beta = 0;
else
  beta = 1/(abs(alpha)*(abs(alpha) + abs(v(1))));
end
v(1) = v(1) + alpha;
x1 = -alpha;