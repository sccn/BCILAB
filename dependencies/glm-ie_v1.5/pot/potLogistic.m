% P = potLogistic(s) - Logistic potential for classification
%
% pot(s) = 1/(1+exp(-s))
%
%   See also POTFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 September 26

function P = potLogistic(s,type,z)

if nargin==1
  q = numel(s); P = zeros(q,4);                                % allocate memory
  P(:,1) = min(0,s(:)) - log(1+exp(-abs(s(:))));          % safe since abs(s)>=0
  p = exp(P(:,1));
  P(:,2) = 1-p;           % 1st derivative of log potential, exp(-s)/(1+exp(-s))
  P(:,3) = p.*(p-1);   % 2nd derivative of log potential, -exp(-s)/(1+exp(-s))^2
  P(:,4) = ones(size(s))/2;
else
  if strcmp(type,'VB')
    P = potLogistic(s);
  elseif strcmp(type,'EP')
    q = numel(s); P = zeros(q,3);                              % allocate memory
    % likLogistic(t) \approx 1/2 + \sum_{i=1}^5 (c_i/2) erf(lam_i/sqrt(2)t)
    lam = sqrt(2)*[0.44 0.41 0.40 0.39 0.36];    % approx coeffs lam_i and c_i
    c = [1.146480988574439e+02; -1.508871030070582e+03; 2.676085036831241e+03;
        -1.356294962039222e+03;  7.543285642111850e+01                      ];
    [lZc,dlZc,d2lZc] = intErf(s*lam, z*(lam.^2));
    P(:,1) = log_expA_x(lZc,c);     % A=lZc, B=dlZc, d=c.*lam', lZ=log(exp(A)*c)
    P(:,2) = expABz_expAx(lZc, c, dlZc, c.*lam');  % ((exp(A).*B)*d)./(exp(A)*c)
      % d2lZ = ((exp(A).*Z)*e)./(exp(A)*c) - dlZ.^2 where e = c.*(lam.^2)'
    P(:,3) = expABz_expAx(lZc, c, dlZc.^2+d2lZc, c.*(lam.^2)') - P(:,2).^2;
    % The scale mixture approximation does not capture the correct asymptotic
    % behavior; we have linear decay instead of quadratic decay as suggested
    % by the scale mixture approximation. By observing that for large values 
    % of -s ln(pot(s)) of potLogistic is linear in s with slope 1, we are
    % able to analytically integrate the tail region; there is no contribution
    % to the second derivative
    val = abs(s)-196/200*z-4;           % empirically determined bound at val==0
    lam = 1./(1+exp(-10*val));                           % interpolation weights
    lZtail = min(z/2-abs(s),-.1);       % apply the same to pot(s) = 1 - pot(-s)
    dlZtail = -sign(s);
    id = s>0; lZtail(id) = log(1-exp(lZtail(id)));        % label and mean agree
    dlZtail(id) = 0;
    P(:,1) = (1-lam).*P(:,1) + lam.* lZtail; % interpolate between scale mixture
    P(:,2) = (1-lam).*P(:,2) + lam.*dlZtail;         % .. and tail approximation    
  else
    error('Unknown type')
  end
end

% Gaussian integral wrt. error function
function [lZ,dlZ,d2lZ] = intErf(mu,s2)
  z = mu./sqrt(1+s2);
  lZ = logphi(z);                                            % log part function
  if nargout>1
    n_p = gauOverCumGauss(z,exp(lZ));
    dlZ = n_p./sqrt(1+s2);                             % 1st derivative wrt mean
    if nargout>2
      d2lZ = -n_p.*(z+n_p)./(1+s2);                    % 2nd derivative wrt mean
    end
  end

% safe implementation of the log of phi(x) = \int_{-\infty}^x N(f|0,1) df
% logphi(z) = log(normcdf(z))
function lp = logphi(z)
  p = (1+erf(z/sqrt(2)))/2;
  lp = zeros(size(z));                                         % allocate memory
  zmin = -6.2; zmax = -5.5;
  ok = z>zmax;                                % safe evaluation for large values
  bd = z<zmin;                                                 % use asymptotics
  ip = ~ok & ~bd;                             % interpolate between both of them
  lam = 1./(1+exp( 25*(1/2-(z(ip)-zmin)/(zmax-zmin)) ));       % interp. weights
  lp( ok) = log( p(ok) );
  % use lower and upper bound acoording to Abramowitz&Stegun 7.1.13 for z<0
  % lower -log(pi)/2 -z.^2/2 -log( sqrt(z.^2/2+2   ) -z/sqrt(2) )
  % upper -log(pi)/2 -z.^2/2 -log( sqrt(z.^2/2+4/pi) -z/sqrt(2) )
  % the lower bound captures the asymptotics
  lp(~ok) = -log(pi)/2 -z(~ok).^2/2 -log( sqrt(z(~ok).^2/2+2)-z(~ok)/sqrt(2) );
  lp( ip) = (1-lam).*lp(ip) + lam.*log( p(ip) );
  
function n_p = gauOverCumGauss(f,p)
  n_p = zeros(size(f));       % safely compute Gaussian over cumulative Gaussian
  ok = f>-5;                            % naive evaluation for large values of f
  n_p(ok) = (exp(-f(ok).^2/2)/sqrt(2*pi)) ./ p(ok); 
  bd = f<-6;                                      % tight upper bound evaluation
  n_p(bd) = sqrt(f(bd).^2/4+1)-f(bd)/2;
  interp = ~ok & ~bd;                % linearly interpolate between both of them
  tmp = f(interp);
  lam = -5-f(interp);
  n_p(interp) = (1-lam).*(exp(-tmp.^2/2)/sqrt(2*pi))./p(interp) + ...
                                                 lam .*(sqrt(tmp.^2/4+1)-tmp/2);

%  computes y = log( exp(A)*x ) in a numerically safe way by subtracting the
%  maximal value in each row to avoid cancelation after taking the exp
function y = log_expA_x(A,x)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  y = log(exp(A-maxA*ones(1,N))*x) + maxA;  % exp(A) = exp(A-max(A))*exp(max(A))
  
%  computes y = ( (exp(A).*B)*z ) ./ ( exp(A)*x ) in a numerically safe way
%  The function is not general in the sense that it yields correct values for
%  all types of inputs. We assume that the values are close together.
function y = expABz_expAx(A,x,B,z)
  N = size(A,2);  maxA = max(A,[],2);      % number of columns, max over columns
  A = A-maxA*ones(1,N);                                 % subtract maximum value
  y = ( (exp(A).*B)*z ) ./ ( exp(A)*x );