% P = potT(s,nu) - Student's t potential
%
% pot(s) = ( 1 + s^2/nu )^(-nu/2-1/2) with nu>0 degrees of freedom
%
%   See also POTFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2011 October 20

function P = potT(s,nu,type,z)

if nargin<=2
  if nargin<2, error('Parameter nu is unspecified.'), end
  q = numel(s); P = zeros(q,4);                           % allocate memory, b=0
  ss = s(:).*s(:);
  P(:,1) = -(nu+1)*log( 1+ss./nu )/2;                            % log potential
  a = s(:) + nu./s(:);
  P(:,2) = -(nu+1)./a;                         % 1st derivative of log potential
  P(:,3) =  (nu+1)*(1-nu./(1e-10+ss))./(a.*a); % 2nd derivative of log potential
else
  if strcmp(type,'VB')
    P = potT(s,nu);
  elseif strcmp(type,'EP')
    error('EP is not supported.')
  else
    error('Unknown type')
  end
end