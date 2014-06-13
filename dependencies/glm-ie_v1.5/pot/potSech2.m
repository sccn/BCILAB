% P = potSech2(s) - Sech-squared potential (Logistic distribution)
%
% pot(s) = 1 / cosh(s)^2
%
%   See also POTFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 September 29

function P = potSech2(s,type,z)

if nargin==1
  q = numel(s); P = zeros(q,4);                           % allocate memory, b=0
  a = exp(-2*abs(s(:))); % always between 0 and 1 and therefore safe to evaluate
  P(:,1) = -2*(abs(s(:)) + log(1+a) - log(2));% log potential, -2*log( cosh(s) )
  b = 2*a./(1+a);
  P(:,2) = -2*sign(s(:)).*(1-b);               % 1st derivative of log potential
  P(:,3) = -4*b./(1+a);                        % 2nd derivative of log potential
else
  if strcmp(type,'VB')
    P = potSech2(s);
  elseif strcmp(type,'EP')
    q = numel(s); P = zeros(q,3);                 % allocate memory, b=0, d2lp=0
    fac = 1e1;          % factor between the widths of the two distributions ...
       % ... from when one considered a delta peak, we use 3 orders of magnitude
    idpot = fac<sqrt(z);                             % Potential is a delta peak
    idgau = sqrt(z)<1/fac;                            % Gaussian is a delta peak
    id = ~idgau & ~idpot;                          % interesting case in between
    % potLogistic(t)   \approx 1/2 + \sum_{i=1}^5 (c_i/2) erf(lam_i/sqrt(2)t)
    % potSech2(t)      \approx \sum_{i=1}^5 c_i potGauss(t|0,rho_i)
    lam = sqrt(2)*[0.44 0.41 0.40 0.39 0.36];  % approx coeffs lam_i, c_i, rho_i
    c   = [1.146480988574439e+02; -1.508871030070582e+03; 2.676085036831241e+03;
          -1.356294962039222e+03;  7.543285642111850e+01                      ];
    rho = sqrt(3)./(pi*lam); o5 = ones(1,5); oq = ones(q,1);
    if any(idpot)
      Q = potGauss(s(idpot)./sqrt(z(idpot))); P(idpot,:) = Q(:,1:3);
      P(idpot,1) = P(idpot,1) - log(2*pi*z(idpot))/2 + log(2);
      P(idpot,2) = P(idpot,2)./sqrt(z(idpot));
      P(idpot,3) = P(idpot,3)./z(idpot);
    end
    if any(idgau)
      Q = potSech2(s(idgau)); P(idgau,:) = Q(:,1:3);
    end
    if any(id)
      mu = s(id)*o5; s2 = z(id)*o5 + oq(id)*rho;              % Gaussian moments
      lZc = -mu.*mu./(2*s2) -log(2*pi*s2)/2 + log(2); 
      dlZc = -mu./s2; d2lZc = -1./s2;
      P(id,1) = log_expA_x(lZc,c);             % A=lZc, B=dlZc, lZ=log(exp(A)*c)
      P(id,2) = expABz_expAx(lZc, c, dlZc, c);     % ((exp(A).*B)*c)./(exp(A)*c)
      P(id,3) = expABz_expAx(lZc, c, dlZc.*dlZc+d2lZc, c) - P(id,2).*P(id,2);

      % the tail asymptotics of potSech2 is the same as for potLaplace
      % which is not covered by the scale mixture approximation, so for
      % extreme values, we approximate potSech2 by a rescaled potLaplace
      tmu = s; tvar = z; crit = 0.596*(abs(tmu)-5.38)-tvar;
      idl = -1<crit & id;                         % if 0<crit, Laplace is better
      if any(idl)                         % close to zero, we use a smooth ..
        lam = 1./(1+exp(-15*crit(idl)));     % .. interpolation with weights lam
        thyp = log(sqrt(6)/pi);
        Q = potLaplace(s(idl)/thyp,type,z(idl));
        P(idl,1) = (1-lam).*P(idl,1) + lam.*Q(:,1);
        P(idl,2) = (1-lam).*P(idl,2) + lam.*Q(:,2)/thyp;
        P(idl,3) = (1-lam).*P(idl,3) + lam.*Q(:,3)/(thyp*thyp);
      end
    end
  else
    error('Unknown type')
  end
end

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