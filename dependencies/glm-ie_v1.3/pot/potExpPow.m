% P = potExpPow(s,al) - Exponential power family potential
%
% pot(s) = exp(-|s|^al) with shape parameter al>0
%
%   See also POTFUNCTIONS.M.
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 September 26

function P = potExpPow(s,al,type,z)

if nargin<=2
  if nargin<2, error('Parameter al is unspecified.'), end
  q = numel(s); P = zeros(q,4);                           % allocate memory, b=0
  P(:,1) = -abs(s(:)).^ al;                                      % log potential
  P(:,2) = -al*sign(s(:)) .*abs(s(:)).^(al-1); % 1st derivative of log potential
  P(:,3) = -al*(al-1)*abs(s(:)).^(al-2);       % 2nd derivative of log potential
else
  if strcmp(type,'VB')
    P = potExpPow(s,al);
  elseif strcmp(type,'EP')
    error('EP is not supported.')
  else
    error('Unknown type')
  end
end