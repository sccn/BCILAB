% Penalised least squares solved by Bregman Splitting. 
% The algorithm is inspired by the code mrics.m from Tom Goldstein 
% http://www.math.ucla.edu/~tagoldst/public_codes/mrics.zip.
%
%   [u,phi] = plsSB(u,X,y,B,opt,lam,pen,varargin)
%
% The method will be fast if two conditions are satisfied
%  (i)  The matrices X'*X and B'*B are diagonal in Fourier space i.e. 
%       the matrix F*A'*A*F' is diagonal for A = {X,B}. This is the
%       case for X, B being equal to matWav, matFFT2line, matFFTNmask, matConv2,
%       matFFTN, matFD2 or scaled vertical concatentations thereof.
%       The function diagFAtAFt.m provides the diagonal matrix.    
%  (ii) The proximity operator pls/prox.m can be evaluated analytically. This is
%       possible for penAbs, penQuad, penZero and penFromPot/potLaplace.
%
% Additional options:
%  opt.
%      SBeta   : coefficient for the constraint term          [default 1]
%                => Data is rescaled such that 1 should work.
%      SBga    : regularization parameter                     [default 0.01/lam]
%      SBinner : number of (inner) loops for the constraint   [default 30]
%                => Sometimes 5-10 iterations are OK.
%      SBncg   : number of (inner) conjugate gradient steps   [default 30]
%                => CG has to be used if the inner update
%                   in u cannot be done analytically
%      SBouter : number of (outer) Bregman iterations         [default 15]
%
%   See also PLSSOLVERS.M.
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 12

function [u,phi] = plsSB(u,X,y,B,opt,lam,pen,varargin)

% optimisation parameters
if isfield(opt,'SBeta')                     % coefficient of the constraint term
  eta = opt.SBeta;
else
  eta = 1;
end
if isfield(opt,'SBga')                                % regularization parameter
  ga = opt.SBga;
else
  ga = 1/(lam*100);
end
if isfield(opt,'SBinner')           % number of (inner) loops for the constraint
  nInner = opt.SBinner;
else
  nInner = 30;
end
if isfield(opt,'SBncg')                    % number of CG steps to solve step a)
  nCG = opt.SBncg;
else
  nCG = 50;
end
if isfield(opt,'SBouter')                 % number of (outer) Bregman iterations
  nOuter = opt.SBouter;
else
  nOuter = 15;
end

% information parameters
if isfield(opt,'output')              % flag saying whether some output is shown
  output = opt.output;
else
  output = false;
end

[q,n] = size(B);
if output, if exist('fflush','builtin'), fflush(stdout); end, end
fac = lam*sqrt(n)/norm(y);  y = fac*y;                      % normalize the data
y0 = y;  s = zeros(q,1);  b = zeros(q,1);                      % allocate memory

% Given that the matrices X and B can be diagonalised in Fourier space i.e. both 
% F*X'*X'*F and F*B'*B*F' are diagonal if F is the Fourier matrix, the
% least-squares problem can be solved by only two ffts and some pointwise 
% operations. Otherwise, we have to run conjugate gradients which is much more
% expensive.
fftdiag = ~isnumeric(X) && ~isnumeric(B);
if fftdiag
  dX = diagFAtAFt(X); dB = diagFAtAFt(B);
  fftdiag = numel(dX)>0 && numel(dB)>0;  
  if fftdiag
    d = dX + eta*dB + ga;
    F = matFFTNmask(true(size(d))); d = d(:);    % orthonormal Fourier transform
  end
else
  if isnumeric(X) && isnumeric(B)
    A = X'*X + eta*B'*B;
  end
end

for o = 1:nOuter
  Xty = X'*y;
  uold = u;
  for i = 1:nInner
    % max_b min_u,s phi(b,u,s) where
    % phi(b,u,s) = ||X*u-y||_2^2 + 2*lam*sum( pen(s) )
    %              + 2*eta*b'*(B*u-s) + eta*||B*u-s||_2^2 + ga*||u-v||_2^2
    % u:   optimisation variable
    % s:   constraint variable
    % b:   Bregman parameters
    % v:   proximity regulariser, here v is the last estimate for u

    % a) update u for b and s fixed
    %    u = argmin 1/2*||X*u-y||_2^2 + eta/2*||B*u-s+b||_2^2 + ga/2*||u-v||_2^2
    %      = A\g,  g = X'*y + eta*B'*(s-b) + ga*v,  A = X'*X + eta*B'*B + ga*I
    g = Xty + eta*B'*(s-b);
    if fftdiag
      if numel(u)>numel(d)
        u = cx2re(F'*((F*re2cx(g+ga*u))./d));                     % make complex
      else
        u = F'*((F*(g+ga*u))./d);
      end
    else
      if isnumeric(X) && isnumeric(B)
        u = A\g;
      else
        u = linsolve_lcg(X,1,B,eta,g,nCG);                 % A = X'*X + eta*B'*B
      end
    end
    
    % b) update s for u and b fixed
    %    s_j <-argmin lam/eta*pen(s_j) + (r_j-s_j)^2/2, r = B*u
    Bu = B*u; r = Bu + b;
    if any(imag(r)>0)                                 % apply proximity operator
      s =    prox(real(r),lam/eta,pen,varargin{:}) ...     % treat real and imag
        + 1i*prox(imag(r),lam/eta,pen,varargin{:});                 % separately
    else
      s = prox(r,lam/eta,pen,varargin{:});
    end
    
    % c) update b for u and s fixed
    %    b <- b + B*u -s
    bold = b;
    b = r - s;
    db = norm(bold-b);
    if db<1e-10, break, end                                          % converged
  end     
  y = y + y0 - X*u;
  du = norm(uold-u)/max(norm(uold),1e-9);
  if output
    phi = feval('phi',u/fac,X,y0/fac,B,lam,pen,varargin{:});
    fprintf('%5i, phi=%4.4e;  du=%1.3e\r', o, phi, du )
    if exist('fflush','builtin'), fflush(stdout); end
  end
  if du<sqrt(eps), break, end                                        % converged
end
u = u/fac;                                                  % undo normalization
if nargout>1, phi = feval('phi',u,X,y0/fac,B,lam,pen,varargin{:}); end
if output, fprintf('\n'); if exist('fflush','builtin'), fflush(stdout); end, end