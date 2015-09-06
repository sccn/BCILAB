function [reg_min,G,reg_param,minG] = gcv(U,s,b,method)
%GCV Plot the GCV function and find its minimum.
%
% [reg_min,G,reg_param] = gcv(U,s,b,method)
% [reg_min,G,reg_param] = gcv(U,sm,b,method)  ,  sm = [sigma,mu]
%
% Plots the GCV-function
%          || A*x - b ||^2
%    G = -------------------
%        (trace(I - A*A_I)^2
% as a function of the regularization parameter reg_param.
% Here, A_I is a matrix which produces the regularized solution.
%
% The following methods are allowed:
%    method = 'Tikh' : Tikhonov regularization   (solid line )
%    method = 'tsvd' : truncated SVD or GSVD     (o markers  )
%    method = 'dsvd' : damped SVD or GSVD        (dotted line)
% If method is not specified, 'Tikh' is default.
%
% If any output arguments are specified, then the minimum of G is
% identified and the corresponding reg. parameter reg_min is returned.

% Per Christian Hansen, IMM, Dec. 16, 2003.

% Reference: G. Wahba, "Spline Models for Observational Data",
% SIAM, 1990.

% Set defaults.
if (nargin==3), method='Tikh'; end  % Default method.
npoints = 200;                      % Number of points on the curve.
smin_ratio = 16*eps;                % Smallest regularization parameter.

% Initialization.
[m,n] = size(U); [p,ps] = size(s);
beta = U'*b; beta2 = norm(b)^2 - norm(beta)^2;
if (ps==2)
  s = s(p:-1:1,1)./s(p:-1:1,2); beta = beta(p:-1:1);
end
if (nargout > 0), find_min = 1; else find_min = 0; end

if (strncmp(method,'Tikh',4) | strncmp(method,'tikh',4))

  % Vector of regularization parameters.
  reg_param = zeros(npoints,1); G = reg_param; s2 = s.^2;
  reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end

  % Intrinsic residual.
  delta0 = 0;
  if (m > n && beta2 > 0), delta0 = beta2; end

  % Vector of GCV-function values.
  for i=1:npoints
    G(i) = gcvfun(reg_param(i),s2,beta(1:p),delta0,m-n);
  end 

  if nargout == 0
      % Plot GCV function.
      loglog(reg_param,G,'-'), xlabel('\lambda'), ylabel('G(\lambda)')
      title('GCV function')
  end

  % Find minimum, if requested.
  if (find_min)
    [minG,minGi] = min(G); % Initial guess.
    reg_min = fminbnd('gcvfun',...
      reg_param(min(minGi+1,npoints)),reg_param(max(minGi-1,1)),...
      optimset('Display','off'),s2,beta(1:p),delta0,m-n); % Minimizer.
    minG = gcvfun(reg_min,s2,beta(1:p),delta0,m-n); % Minimum of GCV function.
    if nargout == 0
        ax = axis;
        HoldState = ishold; hold on;
        loglog(reg_min,minG,'*r',[reg_min,reg_min],[minG/1000,minG],':r')
        title(['GCV function, minimum at \lambda = ',num2str(reg_min)])
        axis(ax)
        if (~HoldState), hold off; end
    end
  end

elseif (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4))
   
  % Vector of GCV-function values.
  rho2(p-1) = abs(beta(p))^2;
  if (m > n & beta2 > 0), rho2(p-1) = rho2(p-1) + beta2; end
  for k=p-2:-1:1, rho2(k) = rho2(k+1) + abs(beta(k+1))^2; end
  G = zeros(p-1,1);
  for k=1:p-1
    G(k) = rho2(k)/(m - k + (n - p))^2;
  end
  reg_param = (1:p-1)';

  % Plot GCV function.
  semilogy(reg_param,G,'o'), xlabel('k'), ylabel('G(k)')
  title('GCV function')

  % Find minimum, if requested.
  if (find_min)
    [minG,reg_min] = min(G);
    ax = axis;
    HoldState = ishold; hold on;
    semilogy(reg_min,minG,'*r',[reg_min,reg_min],[minG/1000,minG],':r')
    title(['GCV function, minimum at k = ',num2str(reg_min)])
    axis(ax);
    if (~HoldState), hold off; end
  end

elseif (strncmp(method,'dsvd',4) | strncmp(method,'dgsv',4))

  % Vector of regularization parameters.
  reg_param = zeros(npoints,1); G = reg_param;
  reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end

  % Intrinsic residual.
  delta0 = 0;
  if (m > n & beta2 > 0), delta0 = beta2; end

  % Vector of GCV-function values.
  for i=1:npoints
    G(i) = gcvfun(reg_param(i),s,beta(1:p),delta0,m-n,1);
  end

  % Plot GCV function.
  loglog(reg_param,G,':'), xlabel('\lambda'), ylabel('G(\lambda)')
  title('GCV function')

  % Find minimum, if requested.
  if (find_min)
    [minG,minGi] = min(G); % Initial guess.
    reg_min = fminbnd('gcvfun',...
      reg_param(min(minGi+1,npoints)),reg_param(max(minGi-1,1)),...
      optimset('Display','off'),s,beta(1:p),delta0,m-n,1); % Minimizer.
    minG = gcvfun(reg_min,s,beta(1:p),delta0,m-n,1); % Minimum of GCV function.
    ax = axis;
    HoldState = ishold; hold on;
    loglog(reg_min,minG,'*r',[reg_min,reg_min],[minG/1000,minG],':r')
    title(['GCV function, minimum at \lambda = ',num2str(reg_min)])
    axis(ax)
    if (~HoldState), hold off; end
  end

elseif (strncmp(method,'mtsv',4) | strncmp(method,'ttls',4))

  error('The MTSVD and TTLS methods are not supported')

else
  error('Illegal method')
end