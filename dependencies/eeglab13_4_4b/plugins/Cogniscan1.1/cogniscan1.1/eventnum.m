% eventnum() - converts voltage levels v into integer event number based
% on predefined voltage levels. It also makes sure that discretized value 
% is consistent for at least two samples. However, this can not be tested 
% for first and last value and so those may not be correct. To fix that
% append the previous and next values to v. Voltage values are assumed to
% be calibrated with the corresponding gain.
%
% Usage:
%   >> n = eventnum(v,mark);
%
% Inputs:
%	v 	- voltage levels
%	mark	- predefined voltage levels
%
% Outputs:
%	n	- discrete levels
%
% Authors: Adam Gerson (adg71@columbia.edu, 2004),
%          with Lucas Parra (parra@ccny.cuny.edu, 2004)
%          and Paul Sajda (ps629@columbia,edu 2004)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Adam Gerson, Lucas Parra and Paul Sajda
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
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

function x=eventnum(x,mark)

% make sure vectors have the right orientation
x=x(:);mark=mark(:)';

% convert into integer labels (not the most efficient but uses less memory)
x = repmat(x,[1 size(mark,2)]);
for i=1:length(mark), x(:,i) = x(:,i)-mark(i); end;
x=abs(x); x=x'; [tmp,x]=min(x);

% detect state transitions
transitions = find(x(1:end-2)~=x(2:end-1))+1;

% transitions may take more than one sample so take the following sample
x(transitions) = x(transitions+1);




