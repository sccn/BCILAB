% BV_CLOSE  Closes a TCP connection to BrainVision Recorder
%     BV_CLOSE(H)
%     In:
%     h : handle to an existing BrainVision connection

% Author: Hal Greenwald, The MITRE Corporation, 29-NOV-2011
   
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function bv_close(h)
if ~h.initialized
    return;
end
disp('Cleaning up connection to BrainVision Recorder');
pnet(h.handle, 'close');
if evalin('base', ['exist(' 'sprintf(''%s'', h.name)' ')'])
    evalin('base',['clear ' sprintf('%s',h.name) ';']);
end