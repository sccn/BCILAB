function model = ssm_cycle()

%SSM_CYCLE Create SSMODEL object for cycle component.
%   model = SSM_CYCLE()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

[Tfun Tgra psi{1}]  = fun_cycle();
[Q Qmmask]          = mat_dupvar(1, false, 2);
[Qfun Qgra psi{2}]  = fun_dupvar(1, false, 2, 'omega tilde');
[psi pmask]         = horzcat(psi{:});
model               = [ssm_null ssmodel('sine cycle', 0, [1 0], ssmat(zeros(2), true(2)), eye(2), ssmat(Q, Qmmask), {'T' 'Q'}, [Tfun Qfun], [Tgra Qgra], psi, pmask)];

