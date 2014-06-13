% [p,dp,d2p] = penVBNorm(s,pot,tau,z,G,theta) - Penalty derived from a 
%                                                           norm/group potential
% 
% pen(s) = -log( pot(tau*r) ) where r = sqrt( G*(s^2 + z) ).
%
% Here, theta are additional parameters for the potential function pot
% which is invoked by the call feval(pot, tau.*r, theta{:}).
% 
% Note that d2p is as large as G'*G and has the same sparsity structure.
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 20

function [p,dp,d2p,b,pi] = penVBNorm(s,pot,tau,z,G,varargin)

if nargin<4, z = 1e-6; end             % default value for smoothing parameter z
r = sqrt( G*(s.*s + z) + eps );                 % stabilised z-smoothed variable
if max(imag(r))>0, error('G and z need to be nonnegative.'), end
q = numel(s); qt = numel(r);                                        % dimensions
[lp,dlp,d2lp,b] = colsplit( feval(pot,tau.*r,varargin{:}) );  % eval. potentials
dlp = tau.*dlp; d2lp = tau.*tau.*d2lp;                   % correct for the scale
if any(b)~=0, error('Only symmetric potentials (b=0) are allowed.'), end

p = -lp;                                                               % penalty
if nargout>1
  v = G'*(dlp./r); dp = -s.*v;                                % first derivative
  if nargout>2
    w = dlp./(r.*r.*r) - d2lp./(r.*r);
    GS = G*sparse(1:q,1:q,s,q,q);
    d2p = GS'*sparse(1:qt,1:qt,w,qt,qt)*GS - sparse(1:q,1:q,v,q,q);  % 2nd deriv
  end
end

if nargout>4, pi = G'*abs(-dlp./r); end                      % return optimal pi

function varargout = colsplit(A)               % split a matrix into its columns
ncol = size(A,2); varargout = cell(1,nargout);
for n=1:nargout, if n>ncol, varargout{n}=[]; else varargout{n}=A(:,n); end, end
