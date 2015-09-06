function model = ssm_airline(s)

%SSM_AIRLINE Create SSMODEL object for airline model.
%   model = SSM_AIRLINE([s])
%       s is the seasonal period which defaults to 12.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, s = 12; end
model   = ssm_sarima(0, 1, 1, 0, 1, 1, s);

