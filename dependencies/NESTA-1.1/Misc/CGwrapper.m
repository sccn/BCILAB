function [y,out2] = CGwrapper(A,b,cg_tol,cg_maxit)
% y = CGwrapper(A,b,cg_tol,cg_maxit)
%   is a simple wrapper to MATLAB's "pcg"
%   but will be less verbose
% Also, calling the function without any input
%   arguments will return the number of function calls
%   (and zero-out the counter)
% i.e.
%   [nA,nCG] = CGwrapper()
%   where nA is the number of times A has been called
%   and nCG is the number of times pcg has been called
%   This will also reset both counters
%
% This version has a special feature: the x0 passed to CG
%   is the result of the previous call.  So this should
%   make it more efficient.
%   (Actually, since NESTA calls inv(AA') twice, on two distinct
%   sequences, we keep track of two distinct sequenes of x0}
%
% NESTA Version 1.1
% See also pcg

% Stephen Becker, 11/15/09

persistent nA nCG
persistent x0
persistent evenOdd
if isempty(nA), nA = 0; end
if isempty(nCG), nCG = 0; end
if isempty(x0), x0 = {[],[]}; end
if isempty(evenOdd), evenOdd=false; end

if nargin > 0
    if nargin < 3, cg_tol = 1e-6; end
    if nargin < 4, cg_maxit = 30; end

    [y,flag,relres,iter] = pcg(A,b,cg_tol,cg_maxit,[],[],x0{ evenOdd + 1 } );

    nA = nA + iter;
    nCG = nCG + 1;
    x0{evenOdd+1} = y;

    evenOdd = ~evenOdd;

else
    y = nA;
    out2 = nCG;
    % and reset
    nA = 0;
    nCG = 0;
    x0={[],[]};
end
