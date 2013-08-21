function [J,alpha,beta,T,history] = vbloreta(varargin)
%[J,varargout] = vbloreta(varargin)
%
% Computes the posterior distribution of the parameters J given some data V. 
% The program solves levels of inference: 1) optimization of parameters J, and
% 2) optimization of hyperparameters alpha and beta. See Trujillo-Barreto
% et. al. (2004) for details.
%
% Ut,s2, and iLV are defined as follows: 
%     Y: Nsensors x time points data matrix
%     K: N x P predictor matrix
%     L: sparse P x P square root of the precision matrix 
%     [U,s,V] = svd( K*inv(L) )
%     iLV = inv(L)*V
%     s2  = s.^2
%
% alpha, beta: hyperparameters
% J: estimated parapeters
% 
%                     P(V|J,alpha)*P(J|beta)
% P(J|V,alpha,beta) = ---------------------- 
%                        P(V|alpha,beta)
% 
% Author: Alejandro Ojeda, SCCN/INC/UCSD, Jan-2013
%         Tim Mullen,      SCCN/INC/UCSD, Feb-2013
%
% References:
%   Trujillo-Barreto, N., Aubert-Vazquez, E., Valdes-Sosa, P.A., 2004.
%     Bayesian model averaging in EEG/MEG imaging. NeuroImage 21, 1300???1319

g=arg_define([0 Inf],varargin,...
    arg_norep({'Ut'},[],[],'U'' from SVD'), ...
    arg_norep({'Y','data'},[],[],'Data matrix'), ...
    arg_norep({'s2'},[],[],'S.^2 from SVD'), ...
    arg_norep({'iLV'},[],[],'inv(L)*V'), ...
    arg_norep({'LtV'},[],[],'L''*V'), ...
    arg({'alpha'},[],[],'Alpha hyperparameter. Noise variance'), ...
    arg({'beta'},[],[],'Beta hyperparameter. Spatial variance'), ...
    arg({'maxTol','MaxTolerance'},1e-3,[0 Inf],'Tolerance for hyperparameter update loop','cat','Loreta Options'), ...
    arg({'maxIter','MaxIterations'},100,[1 Inf],'Maximum iterations for hyperparameter update loop','cat','Loreta Options'), ...
    arg({'gridSize','GridSize'},100,[1 Inf],'Lambda grid size.'), ...
    arg({'history','TrackHistory'},false,[],'Track history for hyperparameters'), ...
    arg({'verbose','VerboseOutput'},false,[],'Verbosity','cat','Loreta Options'), ...
    arg({'initNoiseFactor','InitialNoiseFactor'},0.001,[0 Inf],'Fraction of noise level. Used for initializing alpha parameter','cat','Loreta Options'), ...
    arg({'standardize','Standardize'},true,[],'sLORETA. Scale CSD by estimate of noise variance.') ...
    );
        
    
% copy some variables for convenience
Ut      = g.Ut;
Y       = g.Y;
s2      = g.s2;
iLV     = g.iLV;
alpha   = g.alpha;
beta    = g.beta;

% get variable sizes
nch     = length(s2);
nsrc    = size(iLV,1);
ntp     = size(Y,2);
s       = s2.^(0.5);

if g.history
    [history.alpha, ...
     history.beta,  ...
     history.err,   ...
     history.aic] = deal(nan(1,g.maxIter));
else
    history = [];
end
    
% Initialize hyperparameters
if isempty(alpha) || isempty(beta)
    UtY = Ut*Y(:,1);
    tol = max([nch nch])*eps(max(s));
    lambda2 = logspace(log10(tol),log10(max(s)),g.gridSize);
    gcv = zeros(g.gridSize,1);
    for k=1:g.gridSize
        d = lambda2(k)./(s2+lambda2(k));
        f = diag(d)*UtY;
        gcv(k) = dot(f,f,1)/sum(d)^2;
    end
    loc = getMinima(gcv);
    if isempty(loc), loc = 1;end
    loc = loc(end);
    lambda2 = lambda2(loc);
     
    alpha = g.initNoiseFactor*(Y(:)'*Y(:))/nch;
    beta = alpha*lambda2;
end
err = inf;
for it=1:g.maxIter
    
    if err < g.maxTol, break;end   
        
%     % computing statistics
%     H         = Ut'*diag(alpha*s2./(alpha*s2+beta))*Ut;    
%     SSE       = sum( (Y - H*Y).^2 ,2);
%     sigma2    = mean(SSE);    
%     sqrDiag   = sqrt(alpha*s2./(alpha*s2+beta));
%     sqrSj     = bsxfun(@times, LtV, sqrDiag');
%     Sj        = sum(sqrSj.^2,2);
%     Sj        = sqrt(Sj);
%     aic       = -2*log(sigma2) + 2*nsrc;
%     alpha_old = alpha;
%     beta_old  = beta;
    
    % computing statistics
    H = Ut'*diag(alpha.*s2./(alpha.*s2+beta))*Ut;    
    SSE = mean( (Y - H*Y).^2 ,2);
    sigma2 = mean(SSE);    
    kSjk = diag(H);
        
    q = diag(alpha.*s./(alpha.*s2+beta))*Ut;%*Y;
    Sj = mean((iLV*q).^2,2);
    aic = -2*log(sigma2) + 2*nsrc;
    alpha_old = alpha;
    beta_old = beta;
        
    if g.history
        history.alpha(it) = alpha;
        history.beta(it)  = beta;
        history.err(it)   = err;
        history.aic(it)   = aic;
    end
    if g.verbose
        disp([num2str(it) ' => alpha: ' num2str(alpha) '  beta: ' num2str(beta) ' sse: ' num2str(mean(SSE)) ' hyrp. error: ' num2str(err) ' aic: ' num2str(aic)]);
    end
    
    % updating hyperparameters
%     alpha = updateAlpha(SSE,Sj,ntp);
%     beta  = updateBeta(Sj,nsrc);
    
    % updating hyperparameters
    alpha = updateAlpha(SSE,kSjk,ntp);
    beta  = updateBeta(Sj,nsrc);
    
    err = 0.5*abs(sum(alpha_old-alpha)) + 0.5*abs(beta_old-beta); 
%     err = 0.5*abs(alpha_old-alpha) + 0.5*abs(beta_old-beta);
    
end

if it == g.maxIter, fprintf('vbloreta: Maximum iteration reached. Failed to converge.\n');end

% parameter estimation
T = iLV*diag(alpha.*s./(alpha.*s2+beta))*Ut;   % inverse operator
J = T*Y;                                       % current density


% standardized Loreta
if g.standardize
%     E       = sum(Y-H*Y,2);
%     sigma   = E'*E/(nch-trace(H));
%     dT      = 1./sqrt(dot(T,T,2));
%     S       = 1./sigma*dT;
%     J       = bsxfun(@times,J,S);
    
    E = sum(Y-H*Y,2);
    sigma = E'*E/(nch-trace(H));
    dT = 1./sqrt(dot(T,T,2));
    S = 1./sigma*dT;
    S = S./std(eps+S);
    T = bsxfun(@times,T,S);%sqrt(nsrc)*
    J = bsxfun(@times,J,S);%sqrt(nsrc)*
end


if g.history
    history.alpha(it:end) = [];
    history.beta(it:end)  = [];
    history.err(it:end)   = [];
    history.aic(it:end)   = []; 
end
 
end


%---
function indmin = getMinima(x)
fminor = diff(x)>=0;
fminor = ~fminor(1:end-1, :) & fminor(2:end, :);
fminor = [0; fminor; 0];
indmin = find(fminor);
end

% function alpha = updateAlpha(SSE,r,ntp,bp,cp)
% if nargin < 4, bp = 3;end
% if nargin < 5, cp = 3;end
% b = 1./( (1/bp) + 0.5*sum(SSE) + 0.5*sum(r) );
% c = cp + ntp/2;
% alpha = 1/(b*c);
% end

function alpha = updateAlpha(SSE,kSjk,ntp,bp,cp)
if nargin < 4, bp = 3;end
if nargin < 5, cp = 3;end
b = 1./( (1/bp) + 0.5*SSE + 0.5*kSjk );
c = cp + ntp/2;
alpha = median(1./(b*c));
end


% function beta = updateBeta(r,nsrc,bp,cp)
% if nargin < 3, bp = 3;end
% if nargin < 4, cp = 3;end
% b = 1./( 1/bp + sum(r) );
% c = cp+nsrc/2;
% beta = 1/(b*c);
% end

function beta = updateBeta(r,nsrc,bp,cp)
if nargin < 3, bp = 3;end
if nargin < 4, cp = 3;end
b = 1./( 1/bp + sum(r) );
c = cp+nsrc/2;
beta = 1/(b*c);
end
