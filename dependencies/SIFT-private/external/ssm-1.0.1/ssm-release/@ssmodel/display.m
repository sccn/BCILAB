function display(A)

%@SSMODEL/DISPLAY Command window display of state space model.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if isequal(get(0, 'FormatSpacing'), 'compact'), fprintf(1, '%s =\n', inputname(1));
else fprintf(1, '\n%s =\n\n', inputname(1)); end
disp(A);

