% P = potGauss(s) - Gaussian potential
%
% pot(s) = exp(-s^2/2)
%
%   See also POTFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 September 26

function P = potGauss(s,type,z)

if nargin==1
  q = numel(s); P = zeros(q,4);                           % allocate memory, b=0
  P(:,1) = -s(:).*s(:)/2;                                        % log potential
  P(:,2) = -s(:);                              % 1st derivative of log potential
  P(:,3) = -1;                                 % 2nd derivative of log potential
else
  if strcmp(type,'VB')
    P = potGauss(s);
  elseif strcmp(type,'EP')
    q = numel(s); P = zeros(q,3);                              % allocate memory
    P(:,1) = -s(:).*s(:)./(1+z(:))/2 - log(1+z(:))/2;        % log part function
    P(:,2) = -s(:)./(1+z(:));                       % 1st derivative w.r.t. mean
    P(:,3) = -1./(1+z(:));                          % 2nd derivative w.r.t. mean
  else
    error('Unknown type')
  end
end