function exp = exp_not(exp)
% Negate the given expression.
% Output = Not(Input)
%
% In:
%   Input : An expression
%
% Out:
%   Output : The given expression, negated
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-22

if ~exp_beginfun('default','fingerprint_check',0,'fingerprint_create',0,'memoize',0) return; end

exp = ~exp;

exp_endfun;