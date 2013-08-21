% [s,f] = PROX(r,lam,pen,varargin)
%
% PROX - Proximity operator minimising the decoupling criterion f(s) using a
%        Newton algorithm
%
%        f(s) = (s-r).^2/2 + lam*pen(s)
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 November 22

function [s,f] = prox(r,lam,pen,varargin)

if ischar(pen)
  penstr = pen;
else
  penstr = func2str(pen);
  if numel(strfind(penstr,'@'))>0
    penstr = penstr(1+min(strfind(penstr,')')):max(strfind(penstr,'('))-1);
  end
end
if numel(varargin)>0 && strcmp('newton',varargin{end})
  penstr = ''; varargin = varargin(1:end-1);
end
if     strcmp(penstr,'penAbs')           % analytical solution whenever possible
  s = max(abs(r)-lam,0).*sign(r);
  if nargout>1, f = (s-r).^2/2 + lam*abs(s); end
elseif strcmp(penstr,'penQuad')
  s = r/(1+lam);
  if nargout>1, f = (s-r).^2/2 + (lam/2)*(s.*s); end
elseif strcmp(penstr,'penZero')
  s = r;
  if nargout>1, f = (s-r).^2/2; end
elseif strcmp(penstr,'penVB') && numel(varargin)==3 ...
                                   && strcmp(varargin{1},'potLaplace')
  tau = lam.*varargin{2}; z = varargin{3};
  % solve the quartic: x.^4 + a.*x.^3 + b.*x.^2 + c.*x + d = 0
  a = 2*tau; b = tau.*tau - z - r.*r; c = -2*tau.*z; d = -z.*tau.*tau;
  X = solve_quartic(a,b,c,d);
  % x = sqrt(s*s+z) needs to be real and larger than sqrt(z)
  x = X( abs(imag(X))<1e-10 & real(X) > sqrt(z));
  s = sign(r).*sqrt(x.*x-z);
  if nargout>1, f = (s-r).^2/2 + tau*sqrt(s.*s + z); end
else
  args = {lam,pen,varargin{:}};                               % arguments to fun
  nit = 15; nline = 15;                          % #Newton loops, #line searches
  s = newton(r/2,r,nit,nline,@fun,args{:});
  if nargout>1, f = fun(s,r,args{:}); end
end

function varargout = fun(s,r,lam,pen,varargin)   % f(s) = (s-r)^2/2 + lam*pen(s)
  varargout = cell(nargout, 1);  % allocate the right number of output arguments
  [varargout{:}] = feval(pen,s,varargin{:});
  varargout{1} = lam*varargout{1} + (s-r).*(s-r)/2;
  if nargout>1
    varargout{2} = lam*varargout{2} + (s-r);
    if nargout>2
      varargout{3} = lam*varargout{3} + 1;
    end
  end

function s = newton(s,r,nit,nline,fun,varargin)
  q = numel(s);
  ncv = true(q,1);                                          % not converged flag
  i = 0;                                                             % init loop
  th = (r.*r)/1e8;                                       % convergence threshold
  while i<nit && any(ncv)                                     % Newton algorithm
    [f0,df,d2f] = fun(s(ncv),r(ncv),varargin{:});
    d = -df./d2f;                            % Newton direction of not converged
    f = fun(s(ncv)+d,r(ncv),varargin{:});
    oi = f>f0;                    % indices for which the objective did increase
    oiq = false(q,1); oiq(ncv) = oi;                   % oi expanded to length q 
    j = 0; 
    al = ones(q,1);                                                  % stepsizes
    while j<nline && any(oi)
      al(oiq) = al(oiq)/2;                            % reduce step size by half
      f = fun(s(ncv)+al(ncv).*d,r(ncv),varargin{:});
      oi = f>f0; oiq = false(q,1); oiq(ncv) = oi;
      j = j+1;
    end  
    s(ncv) = s(ncv) + al(ncv).*d;          % final Newton step for not converged
    f = fun(s(ncv),r(ncv),varargin{:});
    cv = f0-f<th(ncv) | oi;                                  % converged indices  
    ncv( ncv ) = ~cv;  
    i = i+1;
  end