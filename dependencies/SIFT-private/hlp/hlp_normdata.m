function [normdata] = hlp_normdata(data,method)
%
% Normalize input data to zero mean and unit variance
%
% Input:
%
%   data:       [channels x points x trials] data to normalize
%   method:     (string)    normalization method
%               'ensemble':     subtract the ensemble mean and divide by
%                               ensemble standard deviation
%               'time':         subtract temporal mean and divide by stdev
%
% Out:
%
%   normdata:   [channels x points x trials] matrix of normalized data
%
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
% 
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


[nchs pnts ntr] = size(data);

% ws = warning;
% warning off
if ischar(method)
    method = {method};
end

for k=1:length(method)
    switch lower(method{k})
        case 'ensemble'
            % pointwise subtract ensemble mean (over trials)
            % divide by ensemble stdev
            if ntr==1
                fprintf('multiple trials not available, ignoring ensemble normalization\n');
                normdata = data;
                return;
            end
            normdata = data-repmat(mean(data,3),[1 1 ntr]);
            normdata = normdata./repmat(std(data,0,3),[1 1 ntr]);
            
            fprintf('NORM\n');
        case 'time'
            % pointwise subtract trial mean
            % divide by trial stdev
            normdata = data-repmat(mean(data,2),[1 pnts 1]);
            normdata = normdata./repmat(std(data,0,2),[1 pnts 1]);
    end
end
% warning(ws);


