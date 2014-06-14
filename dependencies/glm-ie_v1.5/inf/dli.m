%   [m,ga,b,z,zu,nlZ,Q,T] = DLI(X,y,s2,B,t,pot,tau,opts,G)
%
% DLI - Double Loop Inference - Bayesian Inference using a Double Loop algorithm
%       The variational criterion whose stationary point is found has the form
%         phi(ga,b) = ln|A| + h(ga,b) + \min_u R(u,ga,b), where
%         A         = X'*X/s2 + B'*diag(1./ga)*B, s=B*u-t,
%         h(ga,b)   = \sum_{j=1}^q h_j(ga_j,b_j)
%         R(u,ga,b) = (1/s2) * ||X*u-y||^2 + s'*diag(1./ga)*s - 2*b'*s.
%
%       The methods differ in the scalar function h_j(ga_j,b_j). We drop the
%       subindex j in the following.
%
%       In the case of the variational bounding technique (VB), we have:
%         h(ga) = max_x -x/ga - 2 \ln pot(sqrt(x)), x = s*s
%
%       Expectation propagation (EP) yields
%         h(ga,b) = -2/eta * (   \ln E[pot(s)^eta]
%                              - \ln E[exp(b*s-pi*s^2/2)^eta] )
%         where E[#] = \int Q_(s) # ds, and Q_(s) = N(s|mu_,rho_) is the cavity
%         distribution with mu_  = (mu-eta*b*rho)/(1-eta*pi*rho) and
%                           rho_ =            rho/(1-eta*pi*rho).
%         Here, mu and rho are the marginal moments of the Gaussian 
%         approximation N(u|m,V).
%
% See also INFENGINE.M.
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 31

function [m,ga,b,z,zu,nlZ,Q,T] = dli(X,y,s2,B,t,pot,tau, opts,G)

if nargin<7, error('We need at least 6 input arguments'), end
if nargin==7, opts = 1; end
if nargin<9, G = 1; else if numel(G)==0, G = 1; end, end           % set default

% set default values
if isfield(opts,'outerZinit')                    % initial value for upper bound
  z0 = opts.outerZinit;
else
  z0 = 0.05;
end
if isfield(opts,'outerGainit')         % initial value for variational parameter
  pi = 1./opts.outerGainit;
else
  pi = 1;
end
if isfield(opts,'outerNiter')                  % number of outer loop iterations
  nout = opts.outerNiter;
else
  nout = 10;
end

if isfield(opts,'outerVarOpts')   % extra params for approx variance computation
  varOpts = opts.outerVarOpts;
else
  varOpts = [];
end
if isfield(opts,'outerMethod')         % method for computing marginal variances
  outerMethod = opts.outerMethod;
else
  outerMethod = 'lanczos';                                             % default
end
switch outerMethod
  case 'woodbury'
    if ~(isnumeric(B) && numel(B)==1 && B==1 && isnumeric(X))
      error('dli: woodbury is not possible');
    end
  case {'full','factorial'}
  case 'lanczos'
    if ~isfield(varOpts,'MVM')
      if isfield(opts,'outerMVM')
        varOpts.MVM = opts.outerMVM;
      else
        varOpts.MVM = 50;
      end
    end
  case 'sample'
    if ~isfield(varOpts,'outerNSamples')
      varOpts.NSamples = 10;
    end
    if ~isfield(varOpts,'outerNcg')
      varOpts.Ncg = 20;
    end
  otherwise
    error('dli: unknown outer method');
end

if isfield(opts,'innerExact')     % Solve inner problem approximately or exactly
  innerExact = opts.innerExact;
else
  innerExact = 0;
end
if isfield(opts,'innerMVM')              % number of inner loop MVMs or CG steps
  ninMVM = opts.innerMVM;
else
  ninMVM = 50;
end
if isfield(opts,'innerType')
  innerType = opts.innerType;
else
  innerType = 'VB';
end

% information parameters
if isfield(opts,'outerOutput')        % flag saying whether some output is shown
  output = opts.outerOutput;
else
  output = false;
end
if ~isfield(opts,'innerOutput')       % flag saying whether some output is shown
  opts.innerOutput = false;
end

n = max(size(X,2),size(B,2)); u = zeros(n,1);
q = numel(B*u);
z = z0(:).*ones(q,1); ldA = 0; nlZ = inf(nout,1);
b = G'*feval(pot,zeros(q,1)); b = b(:,4); % ck: inserted G'* here
PI = 2*acos(0);   % this is the usual pi which is used by the natural parameters

if output, tic, end
for i = 1:nout
  if output, fprintf('BEGIN loop %d/%d\n',i,nout), end
  pi_old = pi; b_old = b;
  
  % 1) INNER loop(u,ga,b)
  if output, fprintf('INNER: %s\n',innerType), end
  it = ['inner',innerType(1:2)]; opts.innerType = innerType;
  % phi0 = phi[z=0] = h(pi,ga) + min_u R(u), R=||X*u-y||^2/s2+s'*(pi.*s)-2*b'*s
  [u,pi,b,phi0]=feval(it,u,pi,b,z,X,y,s2,B,t,pot,tau,G,innerExact,ninMVM,opts);
  
  % 2) NEGATIVE LOG PARTITION FUNCTION
  % -2 log(Z_Q) = min_u R(u,b,pi) + log(det(A)) - C1
  C1 = n*log(2*PI) - size(X,1)*log(2*PI*s2);
  nlZ(i) = ( phi0 + ldA - C1 )/2;
  if output, fprintf('NLZ(%s): %1.4e\n',innerType,nlZ(i)), end
  
  % 3) OUTER loop: compute marginal variances z from X, s2, B, pi
  if i<nout || nargout>3
    if output, fprintf('OUTER: %s\n',outerMethod), end
    switch outerMethod
      case 'woodbury'
        args = {X,1/s2,pi};
      case {'full','factorial'}
        args = {X,1/s2,B,pi};
      case 'lanczos'
        args = {X,1/s2,B,pi,varOpts.MVM,output};
      case 'sample'
        args = {X,1/s2,B,pi,varOpts.NSamples,varOpts.Ncg};
    end
    if nargout>6
      [z,zu,ldA,Q,T] = feval(['diaginv_',outerMethod],args{:});
    else
      [z,zu,ldA] = feval(['diaginv_',outerMethod],args{:});
    end
    clear diaginv 
  end

  % 4) FINISH: recompute the mean and check convergence
  dpi = dev(pi_old,pi); db = dev(b_old,b); cvg = dpi+db<1e-6;
  if output
    fprintf('END loop %d/%d: db=%1.4e, dpi=%1.4e [%1.1fs]\n',i,nout,db,dpi,toc)
  end
  if cvg || i==nout % converged or last loop
    if output && cvg, fprintf('  Converged\n'), end
    m = postmean(X,s2,B,t,pi,y,b,innerExact,ninMVM,u,output); ga = 1./pi;
    break
  end
end

nlZ = nlZ(1:i);

% Inner Loop: Variational Bayes: penalised least squares
function [u,pi,b,phi0]=innerVB(u,pi,b,z,X,y,s2,B,t,pot,tau,G,exact,ninMVM,opts)
  if isfield(opts,'innerVBpls')  % method for solving the PLS inner loop problem
    pls = opts.innerVBpls;
  else
    OCT = exist('OCTAVE_VERSION') ~= 0;       % check if we run Matlab or Octave
    if OCT, pls = 'plsCG'; else pls = 'plsLBFGS'; end % LBFGS and Octave is hard
  end
  pls_opts.output = opts.innerOutput;
  if pls_opts.output, fprintf('  Use %s solver.\n',pls), end
  pls_opts.nMVM = ninMVM;
  if isfield(opts,'innerIt')            % number of inner Newton loops for plsTN
    pls_opts.nit = opts.innerIt;
  else
    pls_opts.nit = 15;
  end
  if strcmp(pls,'plsTN'), pls_opts.exactNewt = exact; end
  if G==1                    % PLS, no ga dependence since implicitely minimised
    [u,phi_z] = feval(pls,u,X,y,B,t,pls_opts,s2,'penVB',pot,tau,z);
  else
    [u,phi_z] = feval(pls,u,X,y,B,t,pls_opts,s2,'penVBNorm',pot,tau,z,G);
  end
  phi0 = phi_z - sum(z.*pi);            % phi0 = phi_[z=0] at the optimum (pi,u)
  if G==1
    [junk,junk,junk,b,pi] = penVB(B*u,pot,tau,z);
  else
    s = B*u; r = sqrt(G*(s.*s+z));
    [lp,dlp] = cols( feval(pot,tau.*r) ); dlp = tau.*dlp;
    pi = G'*abs(-dlp./r); b = zeros(size(pi));
  end

% Inner Loop: Expectation Propagation: mean computation + moment matching
function [u,pi,b,phi0]=innerEP(u,pi,b,z,X,y,s2,B,t,pot,tau,G,exact,ninMVM,opts)
  nsweep = 8;                      % maximum number of sweeps with fixed z = rho
  
  directeffect = 1;
  if directeffect, pi0 = pi; end         % take direct effect of pi into account
  if isfield(opts,'innerEPeta')   % fractional pot{Laplace,Gauss}, eta \in (0,1]
    eta = opts.innerEPeta;
    [eta,tau] = fracEPpar(pot,eta,tau,opts.innerOutput);
  else
    eta = 1;
  end 
  
  for jj=1:nsweep
    pi_old = pi; b_old = b;
    if opts.innerOutput
      str = sprintf('  %0*d/%d',floor(log(nsweep)/log(10)+1),jj,nsweep);
      fprintf(str)
    end
    u = postmean(X,s2,B,t,pi,y,b,exact,ninMVM,u,opts.innerOutput); s = B*u-t;
    irho = 1./z; if directeffect, irho = irho + (pi-pi0); end
    pi_n = irho - eta.*pi;  b_n = s.*irho - eta.*b;          % cavity parameters
    P = feval(pot, tau.*b_n./pi_n, 'EP', (tau.*tau)./pi_n);
    
    % compute criterion i.e. EP free energy
    lZh = P(:,1);   % lZh = \int_s N(s|b_n/pi_n,1/pi_n) pot(s)^eta ds
                    % lZ  = \int_s N(s|b_n/pi_n,1/pi_n) exp(b*s-pi*s^2/2)^eta ds
    aa = 1-eta.*pi./irho; bb = s-eta.*b./irho;
    lZ = ( (s.*s - bb.*bb./aa).*irho + log(aa) )/2;    
    h  = -2*sum( (lZh-lZ)./eta );
    minR = norm(X*u-y)^2/s2 +  s'*(pi.*s) - 2*b'*s;
    phi0 = minR + h; phi_z = phi0+sum(z.*pi);
    
    % beta = dlZ = d lZ / d mu, nu = d2lZ = -d^2 lZ / d mu^2
    dlZ = tau.*P(:,2); d2lZ = tau.*tau.*P(:,3);
    pi  =  (1-eta).*pi -                      d2lZ  ./(1+d2lZ./pi_n);
    pi(pi<0) = 0;               % enforce positivity i.e. lower bound pi by zero
    b   =  (1-eta).*b  + ( dlZ - (b_n./pi_n).*d2lZ )./(1+d2lZ./pi_n);
    if opts.innerOutput                               % display some information
      dd = [dev(pi_old,pi), dev(b_old,b)];
      spc = repmat(' ',[1,length(str)]);
      fprintf('%s  phi=%4.4e; Block moment match.\n',spc,phi_z)
      fprintf('%s  dpi=%1.2e, db=%1.2e\n',spc,dd)
    end    
    if dev(pi_old,pi)+dev(b_old,b)<1e-6, break, end                 % converged?
  end

% set parameters for fractional EP
function [eta,tau] = fracEPpar(pot,eta,tau,output)
  if ischar(pot), pot_str = pot; else pot_str = func2str(pot); end % make string
  eta = min(1,eta); eta = max(eta,1e-5);     % maxe sure, eta is between 0 and 1
  if strfind(pot_str,'potCat')
    P = feval(pot,[]);                 % which potentials have been concatenated
    q = 0; for i=1:length(P.pots), q = max(q,max(P.ids{i})); end
    e0 = eta; eta = e0*ones(q,1);
    if output, fprintf('  fractional Cat:'), end
    for i=1:length(P.pots)
      poi = P.pots{i}; if ~ischar(poi), poi = func2str(poi); end, id = P.ids{i};
      if     strfind(poi,'potLaplace')
        tau(id) = tau(id).*eta(id);       % Laplace: tau <- tau*eta
        if output, fprintf(', Lap[%d] eta=%1.2f',numel(id),e0), end
      elseif strfind(poi,'potGauss')
        tau(id) = tau(id).*sqrt(eta(id)); % Gauss:   tau <- tau*sqrt(eta)
        if output, fprintf(', Gau[%d] eta=%1.2f',numel(id),e0), end
      else
        eta(id) = 1;
        if output, fprintf(', Oth[%d] eta=1',numel(id)), end
      end
    end
    if output, fprintf('\n'), end
  else
    if     strfind(pot_str,'potLaplace')
      tau = tau*eta;       % Laplace: tau <- tau*eta
      if output, fprintf('  fractional Lap, eta=%1.2f\n',eta), end
    elseif strfind(pot_str,'potGauss')
      tau = tau*sqrt(eta); % Gauss:   tau <- tau*sqrt(eta)
      if output, fprintf('  fractional Gau, eta=%1.2f\n',eta), end
    else
      eta = 1;      % no support for fractional other than Gaussian or Laplacian
      if output, fprintf('  full, eta=%1.2f\n',eta), end
    end
  end

function m = postmean(X,s2,B,t,pi,y,b,exact,nMVM,u,output)      % posterior mean
  if nargin<11, output = false; end
  woodbury = isnumeric(B) && numel(B)==1 && isnumeric(X);
  if woodbury, woodbury = woodbury && B==1; end          % Is Woodbury possible?
  d = [X']*y/s2 + [B']*(b+t.*pi);
  if exact
    if woodbury
      if output, fprintf('  Compute exact mean using Woodbury\n'), end
      m = linsolve_woodbury(X,1/s2,pi,d);                     % Woodbury formula
    else
      if output, fprintf('  Compute exact mean using dense matrices\n'), end
      m = linsolve_full(X,1/s2,B,pi,d);                    % use the full matrix
    end
  else
    if output, fprintf('  Approximate mean using LCG\n'), end
    m = linsolve_lcg(X,1/s2,B,pi,d,nMVM,u);
  end

function d = dev(a,b)                                       % relative deviation
  d = norm(a-b)/max([1,norm(a),norm(b)]);

function varargout = cols(A)                   % split a matrix into its columns
ncol = size(A,2); varargout = cell(1,nargout);
for n=1:nargout, if n>ncol, varargout{n}=[]; else varargout{n}=A(:,n); end, end
