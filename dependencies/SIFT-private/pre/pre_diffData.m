
function [data g] = pre_diffData(varargin)
%
% Apply a difference filter to data. Differencing is a standard operation
% to improve stationarity of a time series. A first-order difference filter
% for input X is given by Y(t) = X(t) - X(t-1). This operation can be
% applied repeatedly to obtain an nth-order difference filter [1]
%
% Inputs:
%
%   EEG:        EEG data structure
%
% Optional:     <'Name',Value> pairs
%
%     VerbosityLevel:    Verbosity level. 0 = no output, 1 = text, 2 = graphical  
%                        Possible values: 0,1,2                                   
%                        Default value  : 1                                       
%                        Input Data Type: real number (double)                    
% 
%     DifferencingOrder: Differencing order                                       
%                        Number of times to difference data                       
%                        Input Range  : [0  10]                                   
%                        Default value: 1                                         
%                        Input Data Type: real number (double)  
% Outputs:
%
%   EEG:        processed EEG structure
%   g:          argument specification structure
%
% See Also: pop_pre_prepData(), pre_prepData(), diff()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1 
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
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



g = arg_define([0 1], varargin, ...
        arg_norep('data',mandatory), ...
        arg({'verb','VerbosityLevel'},int32(1),{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
        arg({'difforder','DifferencingOrder'},1,[0 10],'Differencing order. Number of times to difference data')); 

% commit data to workspace
data = g.data;
g = rmfield(g,'data');
    
[nchs npnt ntr] = size(data);

if g.difforder==0
    return;
end

if g.verb
    disp(['Differencing ' num2str(g.difforder) ' times...']); 
end

if g.verb==2
    multiWaitbar('Differencing','Reset','Color',hlp_getNextUniqueColor);
end
    
datmp = zeros(nchs,npnt-g.difforder,ntr);

% difference each trial
for tr=1:ntr
    if g.verb==2
        multiWaitbar('Differencing',tr/ntr);
    end
    datmp(:,:,tr) = diff(data(:,:,tr),g.difforder,2);
end

% differencing reduces the number of datapoints by difforder
% so here we insert difforder random samples at the beginning 
% of each time-series
noiseVar = 0.01*var(datmp(:));
data      = cat(2,noiseVar*randn(nchs,g.difforder,ntr),datmp);

if g.verb==2
    multiWaitbar('Differencing','Close');
end

