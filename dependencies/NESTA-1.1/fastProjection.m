function [x,k,l] = fastProjection( U, S, V, y, b, epsilon, lambda0, DISP )
% [x,niter,lambda] = fastProjection(U, S, V, y, b, epsilon, [lambda0], [DISP] )
%
% minimizes || x - y ||
%   such that || Ax - b || <= epsilon
%
% where USV' = A (i.e the SVD of A)
%
% The optional input "lambda0" is a guess for the Lagrange parameter
%
% Warning: for speed, does not calculate A(y) to see if x = y is feasible
%
% NESTA Version 1.1
%   See also Core_Nesterov

% Written by Stephen Becker, September 2009, srbecker@caltech.edu

DEBUG = true;
if nargin < 8
    DISP = false;
end
% -- Parameters for Newton's method --
MAXIT = 70;
TOL = 1e-8 * epsilon;
% TOL = max( TOL, 10*eps );  % we can't do better than machine precision

m = size(U,1);
n = size(V,1);
mn = min([m,n]);
if numel(S) > mn^2, S = diag(diag(S)); end  % S should be a small square matrix
r = size(S);
if size(U,2) > r, U = U(:,1:r); end
if size(V,2) > r, V = V(:,1:r); end

s = diag(S);
s2 = s.^2;

% What we want to do:
%   b = b - A*y;
%   bb = U'*b;

% if A doesn't have full row rank, then b may not be in the range
if size(U,1) > size(U,2)
    bRange = U*(U'*b);
    bNull = b - bRange;
    epsilon = sqrt( epsilon^2 - norm(bNull)^2 );
end
b = U'*b - S*(V'*y);  % parenthesis is very important!  This is expensive.
    
% b2 = b.^2;
b2 = abs(b).^2;  % for complex data
bs2 = b2.*s2;
epsilon2 = epsilon^2;

% The following routine need to be fast
% For efficiency (at cost of transparency), we are writing the calculations
% in a way that minimize number of operations.  The functions "f"
% and "fp" represent f and its derivative.

% f = @(lambda) sum( b2 .*(1-lambda*s2).^(-2) ) - epsilon^2;
% fp = @(lambda) 2*sum( bs2 .*(1-lambda*s2).^(-3) );
if nargin < 7, lambda0 = 0; end
l = lambda0; oldff = 0;
one = ones(m,1);
alpha = 1;      % take full Newton steps
for k = 1:MAXIT
    % make f(l) and fp(l) as efficient as possible:
    ls = one./(one-l*s2);
    ls2 = ls.^2;
    ls3 = ls2.*ls;
    ff = b2.'*ls2; % should be .', not ', even for complex data
    ff = ff - epsilon2;
    fpl = 2*( bs2.'*ls3 );  % should be .', not ', even for complex data
%     ff = f(l);    % this is a little slower
%     fpl = fp(l);  % this is a little slower
    d = -ff/fpl;
    if DISP, fprintf('%2d, lambda is %5.2f, f(lambda) is %.2e, f''(lambda) is %.2e\n',...
            k,l,ff,fpl ); end
    if abs(ff) < TOL, break; end        % stopping criteria
    l_old = l;
    if k>2 && ( abs(ff) > 10*abs(oldff+100) ) %|| abs(d) > 1e13 )
        l = 0; alpha = 1/2;  
%         oldff = f(0);
        oldff = sum(b2); oldff = oldff - epsilon2;
        if DISP, disp('restarting'); end
    else
        if alpha < 1, alpha = (alpha+1)/2; end
        l = l + alpha*d;
        oldff = ff;
        if l > 0
            l = 0;  % shouldn't be positive
            oldff = sum(b2);  oldff = oldff - epsilon2;
        end
    end
    if l_old == l && l == 0
        if DISP, disp('Making no progress; x = y is probably feasible'); end
        break;
    end
end
% if k == MAXIT && DEBUG, disp('maxed out iterations'); end
if l < 0
    xhat = -l*s.*b./( 1 - l*s2 );
    x = V*xhat + y;
else
    % y is already feasible, so no need to project
    l = 0;
    x = y;
end