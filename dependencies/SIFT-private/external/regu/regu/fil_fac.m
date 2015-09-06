function f = fil_fac(s,reg_param,method,s1,V1)
%FIL_FAC Filter factors for some regularization methods.
%
% f = fil_fac(s,reg_param,method)
% f = fil_fac(sm,reg_param,method)  ,  sm = [sigma,mu]
% f = fil_fac(s,k,'ttls',s1,V1)
%
% Computes all the filter factors corresponding to the
% singular values in s and the regularization parameter
% reg_param, for the following methods:
%    method = 'dsvd' : damped SVD or GSVD
%    method = 'tsvd' : truncated SVD or GSVD
%    method = 'Tikh' : Tikhonov regularization
%    method = 'ttls' : truncated TLS.
% If sm = [sigma,mu] is specified, then the filter factors
% for the corresponding generalized methods are computed.
%
% If method = 'ttls' then the singular values s1 and the
% right singular matrix V1 of [A,b] must also be supplied.
%
% If method is not specified, 'Tikh' is default.

% Per Christian Hansen, IMM, 12/29/97.

% Initialization.
[p,ps] = size(s); lr = length(reg_param);
if (nargin==2), method = 'Tikh'; end
f = zeros(p,lr);

% Check input data.
if (min(reg_param) <= 0)
  error('Regularization parameter must be positive')
end
if ((strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4) | ...
     strncmp(method,'ttls',4)) & max(reg_param) > p)
  error('Truncation parameter too large')
end

% Compute the filter factors.
for j=1:lr
  if (strncmp(method,'cg',2) | strncmp(method,'nu',2) | strncmp(method,'ls',2))
    error('Filter factors for iterative methods are not supported')
  elseif (strncmp(method,'dsvd',4) | strncmp(method,'dgsv',4))
    if (ps==1)
      f(:,j) = s./(s + reg_param(j));
    else
      f(:,j) = s(:,1)./(s(:,1) + reg_param(j)*s(:,2));
    end
  elseif (strncmp(method,'Tikh',4) | strncmp(method,'tikh',4))
    if (ps==1)
      f(:,j) = (s.^2)./(s.^2 + reg_param(j)^2);
    else
      f(:,j) = (s(:,1).^2)./(s(:,1).^2 + reg_param(j)^2*s(:,2).^2);
    end
  elseif (strncmp(method,'tsvd',4) | strncmp(method,'tgsv',4))
    if (ps==1)
      f(:,j) = [ones(reg_param(j),1);zeros(p-reg_param(j),1)];
    else
      f(:,j) = [zeros(p-reg_param(j),1);ones(reg_param(j),1)];
    end
  elseif (strncmp(method,'ttls',4))
    if (nargin==5)
      coef = ((V1(p+1,:).^2)')/norm(V1(p+1,reg_param(j)+1:p+1))^2;
      for i=1:p
        k = reg_param(j);
        f(i,j) = s(i)^2*...
          sum( coef(1:k)./(s1(1:k)+s(i))./(s1(1:k)-s(i)) );
        if (f(i,j) < 0), f(i,j) = eps; end
        if (i > 1)
          if (f(i-1,j) <= eps & f(i,j) > f(i-1,j)), f(i,j) = f(i-1,j); end
        end
      end
    else
      error('The SVD of [A,b] must be supplied')
    end
  elseif (strncmp(method,'mtsv',4))
    error('Filter factors for MTSVD are not supported')
  else
    error('Illegal method')
  end
end