function connmean = stat_getDistribMean(PConn)
% Return the mean of the distribution of an estimator
%
% Inputs:
% 
%       PConn:     connectivity distribution returned from stat_bootstrap()
%
% Outputs:
%
%       connmean:  the mean of the distribution
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

connmethods = hlp_getConnMethodNames(PConn(1));
% otherfields = setdiff_bc(fieldnames(PConn(1)),connmethods);

for cnd=1:length(PConn)
            
    for m=1:length(connmethods)
        % compute mean of distribution across last dimension
        connmean(cnd).(connmethods{m}) = mean(PConn(cnd).(connmethods{m}),ndims(PConn(cnd).(connmethods{m})));   
    end
       
    extrafields = setdiff_bc(fieldnames(PConn(cnd)),connmethods);
    for i=1:length(extrafields)
        connmean(cnd).(extrafields{i}) = PConn(cnd).(extrafields{i});
    end 
end
