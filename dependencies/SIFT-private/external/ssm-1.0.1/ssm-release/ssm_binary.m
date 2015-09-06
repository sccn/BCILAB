function model = ssm_binary()

%SSM_BINARY Create SSMODEL object for binary distribution error model.
%   model = SSM_BINARY()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

model   = ssmodel('binary error', dist_binary, zeros(1, 0), [], [], []);

