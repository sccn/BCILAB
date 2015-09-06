

function varstr = hlp_variableize(strname)
% returns a 'variableized' version of strname with all leading numbers and
% special characters removed
% strname can be a string or a cell array of strings
% NOTE: if strname is a cell array of length 1 then varstr is a string
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if isempty(strname)
    varstr = [];
    return;
end

if ischar(strname)
    strname = {strname};
end

for i=1:length(strname)
    goodchars = ismember_bc(strname{i},num2str((0:9)')) | isletter(strname{i});
    tmp = strname{i}(goodchars);
    if ismember_bc(tmp(1),num2str((0:9)')), tmp(1)=[]; end
    varstr{i} = tmp;
end

if length(varstr)==1
    varstr = varstr{1};
end
