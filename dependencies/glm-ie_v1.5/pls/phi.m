% Penalised least squares objective function and derivatives
%   phi(u) = 1/lambda * ||X*u-y||_2^2 + 2*sum( pen(s) ), s=B*u-t
%
%   [f,df,d2f] = phi(u,X,y,B,t,lam,pen,varargin)  where penaliser is evaluated
%                                               by feval(pen,s,varargin{:})
%
% The return arguments have the following meaning:
%   f   =     phi(u),               [1,1]  scalar
%   df  = d   phi(u) / du, and      [n,1]  gradient vector (same size as u)
%   d2f = d^2 phi(u) / du^2.        [n,n]  Hessian matrix
%
%   See also PLSALGORITHMS.M.
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 30

function [f,df,d2f] = phi(u,X,y,B,t,lam,pen,varargin)

s = B*u-t; r = X*u-y;
if nargout==1
  p = feval(pen,s,varargin{:});
elseif nargout==2
  [p,dp] = feval(pen,s,varargin{:});
else
  [p,dp,d2p] = feval(pen,s,varargin{:});
end
f = norm(r,2)^2/lam + 2*sum(p);
if nargout>1
  df = 2*([X']*r/lam + [B']*dp);
  if nargout>2              % evaluation of the Hessian is dangerous for large n
    if numel(dp)==numel(d2p)
      d2f = 2*(X'*X/lam + B'*diag(d2p)*B);                     % d2p is a vector
    else
      d2f = 2*(X'*X/lam + B'*d2p*B);          % d2p is a scalar or a full matrix
    end
  end
end