function X = hlp_applyscaling(X, si)
% Apply some previously determined scaling structure
% X = hlp_applyscaling(X, ScaleInfo)
% 
% This is just a convenience tool to implement simple data (e.g. feature) scaling operations.
%
% Examples:
%   scaleinfo = hlp_findscaling(data,'whiten')
%   hlp_applyscaling(data,scaleinfo)
% 
% See also:
%   hlp_findscaling
%
%                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-03-28

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if isfield(si,'add') 
    X = X+repmat(si.add,[size(X,1),1]); end
if isfield(si,'mul') 
    X = X.*repmat(si.mul,[size(X,1),1]); end
if isfield(si,'project') 
    X = X*si.project; end
if isfield(si,{'add','mul','project'}) 
    X(~isfinite(X(:))) = 0; end
