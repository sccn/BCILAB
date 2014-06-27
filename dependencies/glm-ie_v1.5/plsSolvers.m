% A PLS (penalised least squares) solver is a program solving the
% minimisation problem
%   phi(u) = 1/lambda * ||X*u-y||_2^2 + 2*sum( pen(s) ), s = B*u-t,
% where pen(s) is a penalty function.
%
% [u,phi] = pls<NAME>(u0,X,y,B,t,opt,lam,pen,varargin)
%
% INPUT
% =====
%    u0  [nx1]  initial vector
%    X   [mxn]  matrix or operator
%    y   [mx1]  vector
%    B   [qxn]  matrix or operator
%    t   [q,1] vector or [1,1] scalar
%    opt.       optimisation parameters
%      nMVM     maximal number of steps = matrix vector multiplications ... 
%                            with A=X'*X/lam+B'*D*B, D diagonal,   [default 100]
%      output   flag saying whether something is printed         [default false]
%               the function will show the current iteration number, the actual
%               function value phi and the current length of the step
%               in the format:                 13, phi=8.6063e-01;  du=2.325e-04
%    lam [1x1]  (positive) weight of the penaliser, can be Inf
%    pen        penaliser function handle or function name
%                                         [p,dp,d2p] = feval(pen,Bu,varargin{:})
%    varargin   additional parameters for pen
%
% OUTPUT
% ======
%    u   [nx1]  optimal solution
%    phi [1x1]  optimal function value
%
% Currently, we have implemented in pls/pls<NAME>.m
%   LBFGS: LIMITED memory BROYDEN-FLETCHER-GOLDFARB-SHANNO
%               quasi Newton or variable metric method
%   TN:    TRUNCATED NEWTON
%               optimisation with CG approximated Newton steps
%   CGBT:  CONJUGATE GRADIENTS with BACKTRACKING line search
%               optimisation using the Armijo rule
%   CG:    CONJUGATE GRADIENTS
%               optimisation with CG code minimize.m by Carl E. Rasmussen
%   SB:    SPLIT BREGMAN
%               optimisation using an augmented Lagrangian approach
%   BB:    BARZILAI/BORWEIN
%               stepsize adjusted gradient method without descent guarantee
%
% Examples:
% >> lam=1; plsLBFGS(u,X,y,B,t,opt,lam,'penQuad')
% >> lam=1; tau=2; z=1.3; plsCG(u,X,y,B,t,opt,lam,'penVB','potLaplace',tau,z)
%
%   See also PENFUNCTIONS.M, POTFUNCTIONS.M.
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 30
