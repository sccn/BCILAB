function [J,lambdaOpt,T] = ridgeSVD(Y,Ut, s2,V,nlambda,plotGCV,verb)
%[J,lambdaOpt,T] = ridgeSVD(Y,Ut, s2,V,nlambda,plotGCV)
%
% Estimates a ridge regression model, also know as Tikhonov regularization, 
% or minimum norm with L2 prior (or Loreta in the EEG inverse solution literature). 
% For an implementation of sLORETA model see the function inverseSolutionLoreta.
%
% Y: measurements (Nsensors X 1)
% Ut, s2,V are defined as the SVD decomposition of the standardized lead field matrix
% nlambda: maximum size of the grid for the hyperparameter lambda, default: 100
% plotGCV: plot the GCV curve (true/false), default: false
% Jest: estimated parapeters
% T: estimated inverse operatormaximum size of the grid for the hyperparameter lambda, default: 100
% 
% Jest = argmin(J) ||Y-K*J||^2 + lambda*||L*J||^2 == argmin(J) ||Y-K/L*Jst||^2 + lambda*||I||^2, s.t. J = L/Jst 
% and lambda > 0
%
% This code is based on a previous implementation used in Valdes-Hernandez 
% et al. (2009), written by Alejandro Ojeda and Pedro Valdez-Hernandez at 
% the Cuban Neuroscience Center in 2009.
% 
% Author: Alejandro Ojeda, SCCN/INC/UCSD, Jul-2012
%
% References:
%   Pedro A. Valdés-Hernández, Alejandro Ojeda, Eduardo Martínez-Montes, Agustín
%       Lage-Castellanos, Trinidad Virués-Alba, Lourdes Valdés-Urrutia, Pedro A.
%       Valdes-Sosa, 2009. White matter White matter architecture rather than 
%       cortical surface area correlates with the EEG alpha rhythm. NeuroImage 49
%       (2010) 2328–2339

if nargin < 4, error('Not enough input arguments.'); end
if nargin < 5 || isempty(nlambda) nlambda = 100; end
if nargin < 6 || isempty(plotGCV) plotGCV = false; end
if nargin < 7 || isempty(verb) verb = false; end

n = size(Ut,1);
p = size(V,1);
s = sqrt(s2);
UtY = Ut*Y;

tol = max([n p])*eps(max(s));
lambda = logspace(log10(tol),log10(max(s)),nlambda);
gcv = zeros(nlambda,1);

beta2 = mean(norms(Y).^2 - norms(UtY).^2);
[n,m] = size(Ut);
delta0 = 0;
  if (m > n && beta2 > 0)
      if verb
        fprintf('m>n criterion met\n');
        fprintf('m=%d, n=%d\n',m,n);
      end
      delta0 = beta2; 
  end
for it=1:nlambda
    gcv(it) = gcvfun2(lambda(it),s2,UtY,delta0,m-n);
end

loc = getMinima(gcv);
if verb
    fprintf('GCV min found at loc=%d | lambda=%0.5g\n',loc,lambda(loc)); end
if isempty(loc), 
    fprintf('GCV did not find a minimum.\n');
    loc = length(lambda);
end
loc = loc(end);
lambdaOpt = lambda(loc);

T = V*diag(s./(s2+lambdaOpt.^2))*Ut;
J = T*Y;                            % J = (K'*K+lambda*L'*L)\K'*Y

% J = bsxfun(@minus,J,median(J));
% J = bsxfun(@rdivide,J,(std(J)+eps));

if plotGCV
    figure;
    semilogx(lambda,gcv)
    xlabel('log-lambda');
    ylabel('GCV');
    hold on;
    plot(lambdaOpt,gcv(loc),'rx','linewidth',2)
    grid on;
end

%---
function indmin = getMinima(x)
fminor = diff(x)>=0;
fminor = ~fminor(1:end-1) & fminor(2:end);
fminor = [0; fminor; 0];
indmin = find(fminor);


function G = gcvfun2(lambda,s2,beta,delta0,mn,dsvd)

% Auxiliary routine for gcv.  PCH, IMM, Feb. 24, 2008.

% Note: f = 1 - filter-factors.
if (nargin==5)
   f = (lambda^2)./(s2 + lambda^2);
else
   f = lambda./(s2 + lambda);
end
G = mean((norms(bsxfun(@times,f,beta)).^2 + delta0)/(mn + sum(f))^2);


