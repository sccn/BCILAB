function f = setfunc(f, funcf)

%@SSFUNC/SETFUNC Update state space function.
%   f = SETFUNC(f, funcf)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

[f.f{f.fmask}]  = funcf{:, 1};
[f.df{f.fmask}] = funcf{:, 2};

