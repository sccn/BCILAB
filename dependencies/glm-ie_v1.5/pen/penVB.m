% [p,dp,d2p] = penVB(s,pot,tau,z,theta) - Potential derived penalty
% 
% pen(s) = tau*b*(r-s) - log( pot(tau*r) ) where r = sign(s)*sqrt(s^2 + z).
%
% Here, theta are additional parameters for the potential function pot
% which is invoked by the call feval(pot, tau.*r, theta{:}).
%
% For z->0, we have r->s and hence pen(s)->-log( pot(tau*s) ) => MAP estimation.
%
%   See also PENFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 October 21

function [p,dp,d2p,b,pi] = penVB(s,pot,tau,z,varargin)

if nargin<4, z = 1e-6; end             % default value for smoothing parameter z

sign_s = 2*(s>=0)-1;                                % strict sign mapping 0 to 1 
r = sign_s.*sqrt(s.*s + z);                         % signed z-smoothed variable
[lp,dlp,d2lp,b] = colsplit( feval(pot,tau.*r,varargin{:}) );  % eval. potentials
dlp = tau.*dlp; d2lp = tau.*tau.*d2lp; b = tau.*b;       % correct for the scale

p   = b.*(r-s) - lp;                                                   % penalty
s_r = (abs(s)+eps)./(abs(r)+eps);         % stabilised s./r to cover z=0 and s=0
dp  = (b-dlp).*s_r - b;                                       % first derivative
r_e = r + sign_s*eps;
z_r3  = z./(r_e.*r_e.*r_e);            % stabilised z./r.^3 to cover z=0 and s=0
d2p = (b-dlp).*z_r3 - d2lp.*s_r.*s_r;                        % second derivative

id = z==0 | abs(s./sqrt(z+eps))>1e10;       % correct asymptotics of s -> +/-Inf
p(id) = -lp(id);  dp(id) = -dlp(id);  d2p(id) = -d2lp(id);

if nargout>4, pi = abs( (b-dlp)./(r+sign_s/1.5e8) ); end     % return optimal pi

function varargout = colsplit(A)               % split a matrix into its columns
ncol = size(A,2); varargout = cell(1,nargout);
for n=1:nargout, if n>ncol, varargout{n}=[]; else varargout{n}=A(:,n); end, end