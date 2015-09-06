function base = hlp_getBaselineDistrib(PConn,baseline)
% Return the distribution of an estimator averaged over a baseline
%
% Inputs:
% 
%       PConn:     connectivity distribution returned from stat_bootstrap()
%       baseline:  [min max] baseline interval relative to event (t=0 sec)
%
% Outputs:
%
%       base:      the baseline distribution (expanded out to same dims as 
%                  PConn)
%
% See Also: stat_bootstrap()
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.8
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT
%
% Author: Tim Mullen, 2011, SCCN/INC, UCSD. 
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

connmethods = setdiff_bc(fieldnames(PConn),{'winCenterTimes','erWinCenterTimes','freqs'});

for cnd=1:length(PConn)
    
    for m=1:length(connmethods)
        sz = size(PConn(cnd).(connmethods{m}));
        
        baseidx = getindex(PConn(cnd).erWinCenterTimes,baseline);

        base = (mean(PConn(cnd).(connmethods{m})(:,:,:,baseidx(1):baseidx(2),:),4));   % collapse baseline over timepoints
        base = repmat(base,[1 1 1 sz(4) 1]);
    end
                
end
