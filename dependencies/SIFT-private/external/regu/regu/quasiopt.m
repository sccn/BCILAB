function [reg_min,Q,reg_param] = quasiopt(U,s,b,method)
%QUASIOPT Quasi-optimality criterion for choosing the reg. parameter
%
% [reg_min,Q,reg_param] = quasiopt(U,s,b,method)
% [reg_min,Q,reg_param] = quasiopt(U,sm,b,method)  ,  sm = [sigma,mu]
%
% Plots the quasi-optimality function Q for the following methods:
%    method = 'Tikh' : Tikhonov regularization   (solid line )
%    method = 'tsvd' : truncated SVD or GSVD     (o markers  )
%    method = 'dsvd' : damped SVD or GSVD        (dotted line)
% If no method is specified, 'Tikh' is default.
%
% If any output arguments are specified, then the minimum of Q is
% identified and the corresponding reg. parameter reg_min is returned.

% Per Christian Hansen, IMM, Feb. 21, 2001.

% Set defaults.
npoints = 200;  % Number of points for 'Tikh' and 'dsvd'.
if (nargin==3), method = 'Tikh'; end   % Default method.

% Initialization.
[p,ps] = size(s);
if (ps==2), s = s(p:-1:1,1)./s(p:-1:1,2); U = U(:,p:-1:1); end
xi = (U'*b)./s;
if (nargout > 0), find_min = 1; else find_min = 0; end

% Compute the quasioptimality function Q.
if (strncmp(method,'Tikh',4) | strncmp(method,'tikh',4))

  % Compute a vector of Q-values.
  Q = zeros(npoints,1); reg_param = Q;
  reg_param(npoints) = s(p);
  ratio = (s(1)/s(p))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end
  for i=1:npoints
    Q(i) = quasifun(reg_param(i),s,xi);
  end
 
 % Find the minimum, if requested.
 if (find_min)
   [minQ,minQi] = min(Q); % Initial guess.
   reg_min = fminbnd('quasifun',...
     reg_param(min(minQi+1,npoints)),reg_param(max(minQi-1,1)),...
     optimset('Display','off'),s,xi); % Minimizer.
   minQ = quasifun(reg_min,s,xi); % Minimum of function.
 end

elseif (strncmp(method,'dsvd',4) | strncmp(method,'dgsv',4))

  % Compute a vector of Q-values.
  Q = zeros(npoints,1); reg_param = Q;
  reg_param(npoints,1) = s(p);
  ratio = (s(1)/s(p))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end
  for i=1:npoints
    Q(i) = quasifun(reg_param(i),s,xi,1);
  end

 % Find the minimum, if requested.
 if (find_min)
   [minQ,minQi] = min(Q); % Initial guess.
   reg_min = fminbnd('quasifun',...
     reg_param(min(minQi+1,npoints)),reg_param(max(minQi-1,1)),...
     optimset('Display','off'),s,xi,1); % Minimizer.
   minQ = quasifun(reg_min,s,xi,1); % Minimum of function.
 end

elseif (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4))

  % Compute the quasi-optimality function.
  Q = abs(xi); reg_param = (1:p)';

  % Find the minimum, if requested.
  if (find_min)
    [minQ,minQi] = min(Q); reg_min = reg_param(minQi);
  end

else
  error('Illegal method')
end

% Plot the function.
if (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4))
  semilogy(reg_param,Q,'o'), xlabel('k'), ylabel('Q(k)')
  title('Quasi-optimality function')
  if (find_min)
    ax = axis;
    HoldState = ishold; hold on;
    semilogy([reg_min,reg_min],[minQ,minQ/1000],'--')
    axis(ax);
    if (~HoldState), hold off; end
    title(['Quasi-optimality function, minimum at ',num2str(reg_min)])
  end
else
  if (strncmpi(method,'tikh',4) | ...
      strncmpi(method,'dsvd',4) | strncmpi(method,'dgsv',4))
    loglog(reg_param,Q), xlabel('\lambda'), ylabel('Q(\lambda)')
  else
    loglog(reg_param,Q,':'), xlabel('\lambda'), ylabel('Q(\lambda)')
  end
  title('Quasi-optimality function')
  if (find_min)
    ax = axis;
    HoldState = ishold; hold on;
    loglog([reg_min,reg_min],[minQ,minQ/1000],'--')
    axis(ax)
    if (~HoldState), hold off; end
    title(['Quasi-optimality function, minimum at ',num2str(reg_min)])
  end
end