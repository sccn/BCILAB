% P = potLaplace(s) - Laplace (double exponential) potential
%
% pot(s) = exp(-|s|)
%
%   See also POTFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 03

function P = potLaplace(s,type,z)

if nargin==1
  q = numel(s); P = zeros(q,4);                   % allocate memory, b=0, d2lp=0
  P(:,1) = -abs(s(:));                                           % log potential
  P(:,2) = 1-2*(s(:)>0);                       % 1st derivative of log potential
else
  if strcmp(type,'VB')
    P = potLaplace(s);
  elseif strcmp(type,'EP')
    q = numel(s); P = zeros(q,3);                 % allocate memory, b=0, d2lp=0
    fac = 1e3;          % factor between the widths of the two distributions ...
       % ... from when one considered a delta peak, we use 3 orders of magnitude
    idpot = fac<sqrt(z);                             % Potential is a delta peak
    idgau = sqrt(z)<1/fac;                            % Gaussian is a delta peak
    id = ~idgau & ~idpot;                          % interesting case in between
    if any(idpot)
      Q = potGauss(s(idpot)./sqrt(z(idpot))); P(idpot,:) = Q(:,1:3);
      P(idpot,1) = P(idpot,1) - log(2*pi*z(idpot))/2 + log(2);
      P(idpot,2) = P(idpot,2)./sqrt(z(idpot));
      P(idpot,3) = P(idpot,3)./z(idpot);
    end
    if any(idgau)
      Q = potLaplace(s(idgau)); P(idgau,:) = Q(:,1:3);
    end
    if any(id)
      % substitution to obtain unit variance, zero mean Laplacian
      tmu = s(id)./sqrt(2); tvar = z(id)./2;

      % an implementation based on logphi(t) = log(normcdf(t))
      zp = (tmu+sqrt(2)*tvar)./sqrt(tvar);
      zm = (tmu-sqrt(2)*tvar)./sqrt(tvar);
      ap =  logphi(-zp)+sqrt(2)*tmu;
      am =  logphi( zm)-sqrt(2)*tmu;
      P(id,1) = logsum2exp([ap,am]) + tvar;

      lqp = -zp.^2/2 - log(2*pi)/2 - logphi(-zp);       % log( N(z)/Phi(z) )
      lqm = -zm.^2/2 - log(2*pi)/2 - logphi( zm);
      dap = -exp(lqp-log(z(id))/2) + 1;
      dam =  exp(lqm-log(z(id))/2) - 1;
                        % ( exp(ap).*dap + exp(am).*dam )./( exp(ap) + exp(am) )
      P(id,2) = expABz_expAx([ap,am],[1;1],[dap,dam],[1;1]);
 
      a = 2./sqrt(z(id));
      bp = 1 - (a - zp./z(id)).*exp(lqp);
      bm = 1 - (a + zm./z(id)).*exp(lqm);
      % d2lZ(id) = ( exp(ap).*bp + exp(am).*bm )./( exp(ap) + exp(am) ) ...
      %            - dlZ(id).^2;
      P(id,3) = expABz_expAx([ap,am],[1;1],[bp,bm],[1;1]) - P(id,2).*P(id,2);
    end
  else
    error('Unknown type')
  end
end

% computes y = log( sum(exp(x),2) ) in a numerically safe way by subtracting 
%  the row maximum to avoid cancelation after taking the exp
%  the sum is done along the rows
function [y,x] = logsum2exp(logx)
  N = size(logx,2);
  max_logx = max(logx,[],2);
  % we have all values in the log domain, and want to calculate a sum
  x = exp(logx-max_logx*ones(1,N));
  y = log(sum(x,2)) + max_logx;

function y = expABz_expAx(A,x,B,z)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  A = A-maxA*ones(1,N);                                 % subtract maximum value
  y = ( (exp(A).*B)*z ) ./ ( exp(A)*x ); 

% safe implementation of the log of phi(x) = \int_{-\infty}^x N(f|0,1) df
% logphi(z) = log(normcdf(z))
function lp = logphi(z)
  lp = zeros(size(z));                                         % allocate memory
  zmin = -6.2; zmax = -5.5;
  ok = z>zmax;                                % safe evaluation for large values
  bd = z<zmin;                                                 % use asymptotics
  ip = ~ok & ~bd;                             % interpolate between both of them
  lam = 1./(1+exp( 25*(1/2-(z(ip)-zmin)/(zmax-zmin)) ));       % interp. weights
  lp( ok) = log( (1+erf(z(ok)/sqrt(2)))/2 );
  % use lower and upper bound acoording to Abramowitz&Stegun 7.1.13 for z<0
  % lower -log(pi)/2 -z.^2/2 -log( sqrt(z.^2/2+2   ) -z/sqrt(2) )
  % upper -log(pi)/2 -z.^2/2 -log( sqrt(z.^2/2+4/pi) -z/sqrt(2) )
  % the lower bound captures the asymptotics
  lp(~ok) = -log(pi)/2 -z(~ok).^2/2 -log( sqrt(z(~ok).^2/2+2)-z(~ok)/sqrt(2) );
  lp( ip) = (1-lam).*lp(ip) + lam.*log( (1+erf(z(ip)/sqrt(2)))/2 );