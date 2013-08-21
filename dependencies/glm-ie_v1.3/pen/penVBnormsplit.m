% [p,dp,d2p] = penVBNormsplit(s,pot,tau,z,G,theta) - Penalty derived from a 
%                                                    norm/group potential
% 
% pen(s) = tau*b*(r-s) - log( pot(tau*r) ) where r = sign(G*s)*sqrt( G*(s^2 + z) ).
%
% Here, theta are additional parameters for the potential function pot
% which is invoked by the call feval(pot, tau.*r, theta{:}).
% 
% Note that d2p is as large as G'*G and has the same sparsity structure.
%
% For z->0, we have r->s and hence pen(s)->-log( pot(tau*s) ) => MAP estimation.
%
%   See also PENFUNCTIONS.M.
%   Based on penVB and penVBnorm, (c) Christian Kothe, SCCN, 2011 October 9

function [p,dp,d2p,b,pi] = penVBnormsplit(s,pot,tau,z,G,varargin)

if nargin<4, z = 1e-6; end             % default value for smoothing parameter z

sign_Gs = 2*(G*s>=0)-1;                             % strict sign mapping 0 to 1 
r_abs = sqrt( G*(s.*s + z) + eps );             % stabilised z-smoothed variable
r = sign_Gs.*r_abs;                                        % endowed with a sign
if max(imag(r))>0, error('G and z need to be nonnegative.'), end
q = numel(s); qt = numel(r);                                        % dimensions
[lp,dlp,d2lp,b] = cols( feval(pot,tau.*r,varargin{:}) );   % evaluate potentials

if all(b == 0)
    
    % all symmetric potentials; vanilla penVBnorm
    dlp = tau.*dlp; d2lp = tau.*tau.*d2lp;               % correct for the scale
    p = -lp;                                                           % penalty
    if nargout>1
        v = G'*(dlp./r_abs); dp = -s.*v;                      % first derivative
        if nargout>2
            w = dlp./(r_abs.*r_abs.*r_abs) - d2lp./(r_abs.*r_abs);
            GS = G*sparse(1:q,1:q,s,q,q);
            d2p = GS'*sparse(1:qt,1:qt,w,qt,qt)*GS - sparse(1:q,1:q,v,q,q);
        end
    end
    
    if nargout>4, pi = [G']*abs(-dlp./r_abs); end               % return optimal pi
elseif all(b ~= 0)
    
    % all non-symmetric potentials; make sure that this is not an error
    if nnz(G - diag(diag(G))) > 0
        error('Grouped variables work only with symmetric potentials (b=0).'), end
    
    % otherwise vanilla penVB
    dlp = tau.*dlp; d2lp = tau.*tau.*d2lp; b = tau.*b;   % correct for the scale
    p   = b.*(r-s) - lp;                                               % penalty
    s_r = (abs(s)+eps)./(abs(r)+eps);     % stabilised s./r to cover z=0 and s=0
    dp  = (b-dlp).*s_r - b;                                   % first derivative
    r_e = r + sign_s*eps;
    z_r3  = z./(r_e.*r_e.*r_e);        % stabilised z./r.^3 to cover z=0 and s=0
    d2p = (b-dlp).*z_r3 - d2lp.*s_r.*s_r;                    % second derivative
    
    id = z==0 | abs(s./sqrt(z+eps))>1e10;   % correct asymptotics of s -> +/-Inf
    p(id) = -lp(id);  dp(id) = -dlp(id);  d2p(id) = -d2lp(id);
    
    if nargout>4, pi = abs( (b-dlp)./(r+sign_s/1.5e8) ); end % return optimal pi
    
else
    
    % mix of symmetric and asymmetric potentials
    sy = b==0;                                % places of symm potential outputs
    as = b~=0;                               % places of asymm potential outputs

    boundaries = diff(b ~= 0);          % boundaries between the symmetry blocks
    dlp = tau.*dlp; d2lp = tau.*tau.*d2lp; b = tau.*b;   % correct for the scale
    
    % 1) calc symmetric portion like penVBnorm
    
    % restrict all variables to the symmmetric part
    sy_s = sum(G(sy,:),1) ~= 0;        % places of symm potential inputs (pre-G)
    lp_sy = lp(sy);
    dlp_sy = dlp(sy);
    d2lp_sy = d2lp(sy);
    r_sy = r_abs(sy);
    qt_sy = nnz(sy);
    q_sy = nnz(sy_s);
    s_sy = s(sy_s);        
    G_sy = G; G_sy(~sy,:) = []; G_sy(:,~sy_s) = [];
    b_sy = b(sy);
    
    p_sy = -lp_sy;                                                     % penalty
    if nargout>1
        v_sy = G_sy'*(dlp_sy./r_sy);
        dp_sy = -s_sy.*v_sy;                                  % first derivative
        if nargout>2
            w_sy = dlp_sy./(r_sy.*r_sy.*r_sy) - d2lp_sy./(r_sy.*r_sy);
            GS_sy = G_sy*sparse(1:q_sy,1:q_sy,s_sy,q_sy,q_sy);
            d2p_sy = GS_sy'*sparse(1:qt_sy,1:qt_sy,w_sy,qt_sy,qt_sy)*GS_sy - sparse(1:q_sy,1:q_sy,v_sy,q_sy,q_sy);
        end
    end
    
    if nargout>4
        pi_sy = [G_sy']*abs(-dlp_sy./r_sy); end              % return optimal pi
        
    
    % 2) calc asymmetric portion like penVB
    
    % restrict all variables to the asymmetric part
    as_s = ~sy_s;                     % places of asymm potential inputs (pre-G)
    lp_as = lp(as);
    dlp_as = dlp(as);
    d2lp_as = d2lp(as);
    r_as = r(as);
    qt_as = nnz(as);
    q_as = nnz(as_s);
    s_as = s(as_s);        
    G_as = G; G_as(~as,:) = []; G_as(:,~as_s) = [];
    b_as = b(as);
    sign_Gs_as = sign_Gs(as);
    z_as = z(as_s);
    if nnz(G_as-diag(diag(G_as))) > 0
        error('When restricted to asymmetric potential places, G must be diagonal.'); end
    
    p_as   = b_as.*(r_as-s_as) - lp_as;                                % penalty
    s_r_as = (abs(s_as)+eps)./(abs(r_as)+eps);% stabilised s./r to cover z=0 and s=0
    dp_as  = (b_as-dlp_as).*s_r_as - b_as;                    % first derivative
    r_e_as = r_as + sign_Gs_as*eps;
    z_r3_as  = z_as./(r_e_as.*r_e_as.*r_e_as);% stabilised z./r.^3 to cover z=0 and s=0
    d2p_as = (b_as-dlp_as).*z_r3_as - d2lp_as.*s_r_as.*s_r_as;% second derivative
    
    id_as = z_as==0 | abs(s_as./sqrt(z_as+eps))>1e10;   % correct asymptotics of s -> +/-Inf
    p_as(id_as) = -lp_as(id_as);  dp_as(id_as) = -dlp_as(id_as);  d2p_as(id_as) = -d2lp_as(id_as);
    
    if nargout>4
        pi_as = abs( (b_as-dlp_as)./(r_as+sign_Gs_as/1.5e8) ); end % return optimal pi
    
    if nnz(boundaries) == 1
        if sum(boundaries) > 0
            % first symmetric followed by asymmetric
            p = [p_sy; p_as];
            if nargout > 1
                dp = [dp_sy; dp_as];
                if nargout > 2
                    d2p = blkdiag(d2p_sy,diag(d2p_as));
                    b = [b_sy; b_as];
                    if nargout > 4
                        pi = [pi_sy; pi_as]; end
                end
            end
        else
            % first asymmetric followed by symmetric            
            p = [p_as; p_sy];
            if nargout > 1
                dp = [dp_as; dp_sy];
                if nargout > 2
                    d2p = blkdiag(diag(d2p_as),d2p_sy);
                    b = [b_as; b_sy];
                    if nargout > 4
                        pi = [pi_as; pi_sy]; end                        
                end
            end
        end
    else
        error('Mixed organization of symmetric and asymmetric potentials not yet supported if grouping enabled.');
    end
end

function varargout = cols(A)                   % split a matrix into its columns
ncol = size(A,2); varargout = cell(1,nargout);
for n=1:nargout, if n>ncol, varargout{n}=[]; else varargout{n}=A(:,n); end, end

