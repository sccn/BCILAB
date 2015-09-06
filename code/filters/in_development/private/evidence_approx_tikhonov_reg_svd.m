function [x,alpha,beta,T] = evidence_approx_tikhonov_reg_svd(y,Ut,s2,iLV,alpha,beta,options)
[n,m] = size(iLV);
if nargin < 5, alpha = rand;end
if nargin < 6, beta = rand;end
if nargin < 7, 
    options = struct('maxTol',1e-6,'maxIter',100,'verbose',false);
end
if ~isfield(options,'maxTol'), options.maxTol = 1e-6;end
if ~isfield(options,'maxIter'), options.maxIter = 100;end
if ~isfield(options,'verbose'), options.verbose = false;end
if isempty(alpha), alpha = rand;end
if isempty(beta), beta = rand;end

s = sqrt(s2);
UtY = Ut*y;

useGPU = isa(y,'gpuArray');
if useGPU
    alpha_k = gpuArray.zeros(options.maxIter,1);
    beta_k = gpuArray.zeros(options.maxIter,1);
    gamma_k = gpuArray.zeros(options.maxIter,1);
else
    alpha_k = zeros(options.maxIter,1);
    beta_k = zeros(options.maxIter,1);
    gamma_k = zeros(options.maxIter,1);
end
alpha_k(1) = alpha;
beta_k(1) = beta;
gamma_k(1) = m;
if options.verbose
    fprintf('Iter\tGamma\tAlpha\t\tBeta\n');
end
for k=2:options.maxIter
    norm_x = mean(sum((bsxfun(@times,beta*s./(beta*s2+alpha),UtY)).^2));
    norm_e = mean(sum(bsxfun(@times,alpha./(beta*s2+alpha),UtY).^2));

    gamma = sum(s2.*beta./(alpha+s2.*beta));
    alpha = gamma/norm_x;
    beta = (m-gamma)/norm_e;
    
    alpha_k(k) = alpha;
    beta_k(k) = beta;
    gamma_k(k) = gamma;
    
    if options.verbose
        fprintf('%i\t%5g\t%e\t%e\n',k,gamma,alpha,beta);
    end
    if any([abs([diff(alpha_k(k-1:k)) diff(beta_k(k-1:k))])] < options.maxTol)
        break;
    end
end
[x,T] = tikhonoc_reg(y,Ut,s,iLV,alpha,beta);
end


function [x,T] = tikhonoc_reg(y,Ut,s,V,alpha,beta)
if nargin < 5, alpha = rand;end
if nargin < 6, beta = rand;end
s2 = s.^2;
T = V*bsxfun(@times,beta*s./(beta*s2+alpha),Ut);
x = T*y;
end