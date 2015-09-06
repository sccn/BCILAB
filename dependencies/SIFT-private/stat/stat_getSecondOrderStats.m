function Stats = stat_getSecondOrderStats(PConn,alpha)
% Return the mean, standard deviation, and confidence intervals of a
% surrogate distribution
%
% Inputs:
% 
%       PConn:     connectivity distribution returned from stat_bootstrap()
%       alpha:     significance level for empircal confidence intervals 
%                  (e.g. alpha=0.95 for 95% confidence intervals)
% Outputs:
%
%       Stats:     Structure containing
%                  .mean 
%                  .stdev
%                  .ci
%                  for each connectivity measure
%
% See Also: stat_surrogate(), stat_surrogateStats()
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

connmethods = hlp_getConnMethodNames(PConn);

for cnd=1:length(PConn)
    
    for m=1:length(connmethods)
        
        sz = size(PConn(cnd).(connmethods{m}));

        Stats(cnd).(connmethods{m}).mean  = mean(PConn(cnd).(connmethods{m}),length(sz));
        Stats(cnd).(connmethods{m}).stdev = std(PConn(cnd).(connmethods{m}),0,length(sz));
        Stats(cnd).(connmethods{m}).ci(1,:,:,:,:,:) = ...
            single(prctile(PConn(cnd).(connmethods{m}),(100*alpha)/2,length(sz)));     % lower ci
        Stats(cnd).(connmethods{m}).ci(2,:,:,:,:,:) = ...
            single(prctile(PConn(cnd).(connmethods{m}),100-100*alpha/2,length(sz)));   % upper ci
        
    end
                
end
            