function base = stat_getBaselineDistrib(PConn,baseline,winCenterTimes,noavg)
% Return the distribution of an estimator averaged over a baseline and
% expanded out to the same dimensions as matrices in PConn connectivity
% object
%
% Inputs:
% 
%       PConn:     connectivity distribution returned from stat_bootstrap()
%       baseline:  [min max] baseline interval relative to event (t=0 sec)
%
% Optional
%
%       winCenterTimes: vector of times corresponding to window centers
%                      (same units as baseline)
%       noavg:     if true, the baseline distribution is returned 
% Outputs:
%
%       base:      the baseline mean distribution (expanded out to same dims as 
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

if nargin<4
    noavg = false;
end

if isstruct(PConn)
    % connectivity object
    connmethods = hlp_getConnMethodNames(PConn);

    for cnd=1:length(PConn)

        baseidx = getindex(PConn(cnd).erWinCenterTimes,baseline);

        for m=1:length(connmethods)
            sz = size(PConn(cnd).(connmethods{m}));
            if noavg
                base.(connmethods{m}) = PConn(cnd).(connmethods{m})(:,:,:,baseidx(1):baseidx(2),:);
            else
                basetmp = (mean(PConn(cnd).(connmethods{m})(:,:,:,baseidx(1):baseidx(2),:),4));  % collapse baseline over timepoints
                base.(connmethods{m}) = repmat(basetmp,[1 1 1 sz(4) 1]);  % expands out to appropriate size
            end
        end

    end
    
else
    if nargin<3
        error('SIFT:stat_getBaselineDistrib','You must provide winCenterTimes');
    end
    
    % single matrix
    baseidx = getindex(winCenterTimes,baseline);
    sz = size(PConn);
    if noavg
        base = PConn(:,:,:,baseidx(1):baseidx(2),:);
    else
        basetmp = (mean(PConn(:,:,:,baseidx(1):baseidx(2),:),4));   % collapse baseline over timepoints
        base = repmat(basetmp,[1 1 1 sz(4) 1]);  % expands out to appropriate size
    end
end

    
