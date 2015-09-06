function model = ssm_const()

%SSM_CONST Create SSMODEL object for constant shift.
%   model = SSM_CONST()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

model   = [ssm_null ssmodel('constant', 0, 1, 1, zeros(1, 0), [])];

