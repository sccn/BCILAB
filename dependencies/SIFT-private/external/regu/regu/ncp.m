function [reg_min,dist,reg_param] = ncp(U,s,b,method)
%NCP Plot the NCPs and find the one closest to a straight line.
%
% [reg_min,G,reg_param] = ncp(U,s,b,method)
% [reg_min,G,reg_param] = ncp(U,sm,b,method)  ,  sm = [sigma,mu]
%
% Plots the normalized cumulative priodograms (NCPs) for the residual
% vectors A*x - b.  The following methods are allowed:
%    method = 'Tikh' : Tikhonov regularization
%    method = 'tsvd' : truncated SVD or GSVD
%    method = 'dsvd' : damped SVD or GSVD
% If method is not specified, 'Tikh' is default.
%
% The NCP closest to a straight line is identified and the corresponding
% regularization parameter reg_min is returned.  Moreover, dist holds the
% distances to the straight line, and reg_param are the corresponding
% regularization parameters.

% Per Christian Hansen, IMM, Jan. 4, 2008.

% Reference: P. C. Hansen, M. Kilmer & R. H. Kjeldsen, "Exploiting
% residual information in the parameter choice for discrete ill-posed
% problems", BIT 46 (2006), 41-59.

% Set defaults.
if (nargin==3), method='Tikh'; end  % Default method.
npoints = 200;                      % Number of initial NCPS for Tikhonov.
nNCPs = 20;                         % Number of NCPs shown for Tikhonov.
smin_ratio = 16*eps;                % Smallest regularization parameter.

% Initialization.
m = size(U,1); [p,ps] = size(s);
beta = U'*b;
if (ps==2)
  s = s(p:-1:1,1)./s(p:-1:1,2); beta = beta(p:-1:1);
end

if (strncmp(method,'Tikh',4) | strncmp(method,'tikh',4))

  % Vector of regularization parameters.
  reg_param = zeros(npoints,1);
  reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end

  % Vector of distances to straight line.
  dists = zeros(npoints,1);
  if isreal(beta), q = floor(m/2); else q = m-1; end
  cp = zeros(q,npoints);
  for i=1:npoints
    [dists(i),cp(:,i)] = ncpfun(reg_param(i),s,beta(1:p),U(:,1:p));
  end 

  % Plot selected NCPs.
  stp = round(npoints/nNCPs);
  plot(cp(:,1:stp:npoints)), hold on

  % Find minimum.
  [minG,minGi] = min(dists); % Initial guess.
  reg_min = fminbnd('ncpfun',...
    reg_param(min(minGi+1,npoints)),reg_param(max(minGi-1,1)),...
    optimset('Display','off'),s,beta(1:p),U(:,1:p)); % Minimizer.
  [dist,cp] = ncpfun(reg_min,s,beta(1:p),U(:,1:p));
  plot(cp,'-r','linewidth',3), hold off
  title(['Selected NCPs.  Most white for \lambda = ',num2str(reg_min)])

elseif (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4))
   
  % Matrix of residual vectors.
  R = zeros(m,p-1);
  R(:,p-1) = beta(p)*U(:,p);
  for i=p-1:-1:2
      R(:,i-1) = R(:,i) + beta(i)*U(:,i);
  end
  
  % Compute NCPs and distances.
  if isreal(R), q = floor(m/2); else q = m-1; end
  D = abs(fft(R)).^2; D = D(2:q+1,:);
  v = (1:q)'/q; cp = zeros(q,p-1); dist = zeros(p-1,1);
  for k=1:p-1
    cp(:,k) = cumsum(D(:,k))/sum(D(:,k));
    dist(k) = norm(cp(:,k)-v);
  end

  % Locate minimum and plot.
  [dist_min,reg_min] = min(dist);
  plot(cp), hold on
  plot(1:q,cp(:,reg_min),'-r','linewidth',3), hold off
  title(['Most white for k = ',num2str(reg_min)])
  
  reg_param = (1:p-1)';
  
elseif (strncmp(method,'dsvd',4) | strncmp(method,'dgsv',4))

  % Vector of regularization parameters.
  reg_param = zeros(npoints,1);
  reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = ratio*reg_param(i+1); end

  % Vector of distances to straight line.
  dists = zeros(npoints,1);
  if isreal(beta), q = floor(m/2); else q = m-1; end
  cp = zeros(q,npoints);
  for i=1:npoints
    [dists(i),cp(:,i)] = ncpfun(reg_param(i),s,beta(1:p),U(:,1:p),1);
  end 

  % Plot selected NCPs.
  stp = round(npoints/nNCPs);
  plot(cp(:,1:stp:npoints)), hold on

  % Find minimum, if requested.
  [minG,minGi] = min(dists); % Initial guess.
  reg_min = fminbnd('ncpfun',...
    reg_param(min(minGi+1,npoints)),reg_param(max(minGi-1,1)),...
    optimset('Display','off'),s,beta(1:p),U(:,1:p),1); % Minimizer.
  [dist,cp] = ncpfun(reg_min,s,beta(1:p),U(:,1:p));
  plot(cp,'-r','linewidth',3), hold off
  title(['Selected NCPs.  Most white for \lambda = ',num2str(reg_min)])

elseif (strncmp(method,'mtsv',4) | strncmp(method,'ttls',4))
  error('The MTSVD and TTLS methods are not supported')
else
  error('Illegal method')
end