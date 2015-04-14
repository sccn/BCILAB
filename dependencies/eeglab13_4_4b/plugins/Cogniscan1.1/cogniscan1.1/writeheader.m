% writeheader() - Write CogniScan header file
%
% Usage:
%   >> writeheader(file,D,N,fs,gain)
%
% Inputs:
%   file        - name of CogniScan data file
%   D           - Number of channels
%   N           - Number of samples (i.e. slices)
%   fs          - sampling rate (Hz)
%   gain        - conversion factor to convert data to Volts
%
% Authors: Adam Gerson (adg71@columbia.edu, 2004),
%          with Lucas Parra (parra@ccny.cuny.edu, 2004)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Adam Gerson and Lucas Parra
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

function writeheader(file,D,N,fs,gain)

fid = fopen([file '.hdr'],'w');
fprintf(fid, 'Slice                          ; Mode\n');
fprintf(fid, '%d                             ; Channels\n', D);
fprintf(fid, '%d                             ; Slices\n', N);
fprintf(fid, '%e                             ; Sample period\n', 1/fs);
fprintf(fid, '%e                             ; Conversion Factor\n', gain);
fclose(fid);
