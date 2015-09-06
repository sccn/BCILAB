function Conn2 = hlp_rmbaseline(Conn,baseline,varargin)
% remove the baseline from connectivity matrices
%
% Input: 
%
%   Conn:       Either a 4D connectivity matrix (nchs x nchs x freq x time)
%               or a Connectivity object (structure) with fields
%               winCenterTimes and connectivity methods
%   baseline:   The baseline to remove. This should be [min max] in seconds 
%               relative to winCenterTimes (e.g., [0 0.5] for first 500ms)
%   varargin:   - If Conn is a structure, then this is an (optional) cell
%               array of connectivity measures to remove baseline from
%               (subset of Conn)
%               - If Conn is a matrix, this is a *mandatory* vector of
%               winCenterTimes
%
% Output:
%
%   Conn2:      Connectivity object with mean baseline removed from each
%               connectivity measure
%
% See Also: est_mvarConnectivity()
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
if isstruct(Conn)
    % Conn is connectivity object

    if nargin<3
        connmethods = hlp_getConnMethodNames(Conn);
    elseif iscell(varargin{1})
        connmethods = varargin{1};
    else
        error('incorrect option for ''connmethods''');
    end

    baseidx = getindex(Conn.winCenterTimes,baseline);

    for m=1:length(connmethods)
        Conn2.(connmethods{m}) = rmbase(Conn.(connmethods{m}),baseidx);
    end
    
    % copy supplementary fields into new connectivity structure
    extrafields = setdiff_bc(fieldnames(Conn),hlp_getConnMethodNames(Conn));
    for f=1:length(extrafields)
        Conn2.(extrafields{f}) = Conn.(extrafields{f});
    end
    
else
    % Conn is a 4-dimensional matrix
    if nargin < 3
        error('you must provide an array of winCenterTimes as the third argument');
    else
        winCenterTimes = varargin{1};
    end

    baseidx = getindex(winCenterTimes,baseline);

    Conn2 = rmbase(Conn,baseidx);
end
        
        
function debased = rmbase(data,baseidx)
    debased = data - repmat(mean(data(:,:,:,baseidx(1):baseidx(end)),4), ...
                        [1,1,1,size(data,4)]); 
                    