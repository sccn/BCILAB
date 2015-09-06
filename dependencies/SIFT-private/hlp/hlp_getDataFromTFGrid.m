function Conn = hlp_getDataFromTFGrid(figh,curComponentNames)
% extract the image data from a Time-Frequency Grid and store in Conn object
% Inputs: 
%   figh:               Time-Frequency Grid figure handle
%   curComponentNames:  cell array containing the names of the variables on
%                       the columns/rows of the grid. This must be the same
%                       as the names stored in EEG.CAT.curComponentNames
%                       used at the time of generating the figure
%
% Example: Conn = hlp_getDataFromTFGrid(gcf,EEG.CAT.curComponentNames);
%
% See Also: vis_TimeFreqGrid()
%
%
% Author: Tim Mullen March, 2012, SCCN/INC, UCSD. 
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

causalplots = findobj(figh,'tag','causalplot');


for handle=causalplots(:)'
    % get the from/to labels for this plot
    usrdata = get(handle,'UserData');
    
    % get the nodelabels and look up [from, to] indices
    for k=1:length(usrdata.nodelabels)
        varidx(k) = find(strcmp(curComponentNames,usrdata.nodelabels{k}));
    end
    
    
    connmethod = usrdata.connmethod;
    
    % get handle of image and extract Cdata
    Conn.(connmethod)(varidx(2),varidx(1),:,:) = get(findobj(handle,'type','image'),'Cdata');
    Conn.erWinCenterTimes = usrdata.alltimes;
    Conn.winCenterTimes   = usrdata.alltimes;
    Conn.freqs            = usrdata.allfreqs;
    
end