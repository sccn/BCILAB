function EEG = hlp_blockAverageConnectivity(EEG,numblocks)
    % INPUTS: 
    % 
    %   EEG         - EEG structure containing EEG.CAT.Conn
    %
    % Optional:
    %
    %   numblocks   - number of blocks to average over
    %
    % OUTPUTS: 
    %
    %   EEG     - contained modifed Conn structure 
    
    % average connectivity across trials or blocks
    % trials/blocks are assumed to be concatenated columnwise in last
    % dimension of EEG.CAT.Conn.<connmethod>
    % e.g. EEG.CAT.Conn.<connmethod> is assumed to be [nchs x nchs x nfreqs x numblocks*blocklength]
    % numblocks must exactly devide the number of time points in Conn matrices
    
    % References:
    %
    % [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
    %     Theoretical Handbook and User Manual. Chapter 4
    %     Available at: http://www.sccn.ucsd.edu/wiki/Sift
    %
    % Author: Tim Mullen 2011, SCCN/INC, UCSD.
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

    
    connmethod = hlp_getConnMethodNames(EEG.CAT.Conn);
    
    if nargin < 2
        % default number of blocks is number of trials
        numblocks = size(EEG.CAT.srcdata,3);
    end
    
    [nchs nchs nfreq npoints] = size(EEG.CAT.Conn.(connmethod{1}));
        
    for m =1:length(connmethod)
        C = EEG.CAT.Conn.(connmethod{m});
        C = reshape(C,[nchs nchs nfreq npoints/numblocks numblocks]);
        EEG.CAT.Conn.(connmethod{m}) = mean(C,5);
        
    end
    
    EEG.CAT.Conn.erWinCenterTimes = EEG.CAT.Conn.erWinCenterTimes(1:npoints/numblocks);
    
    