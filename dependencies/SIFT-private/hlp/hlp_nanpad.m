function data = hlp_nanpad(data,morder)
%
% Transform 3D [points x channels x trials] matrix into 2D
% [(points+morder+2)*trials x channels] matrix with morder+2 NaNs between 
% each trial.
%
% Input: 
% 
%   data        [pnts x chs x trials]
%   morder      VAR model order
%
% Output: 2D nan-padded data
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

[pnts nch ntr] = size(data);

nanc=nan;
data = cat(1,data,nanc(ones(morder+2,nch,ntr))); 
data = permute(data,[1 3 2]); % need it in [pnts,trials,chns]
data = reshape(data,(pnts+morder+2)*ntr,nch);