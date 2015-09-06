function figh=hlp_viewGraphicsResource(rscName,hax)
% render a SIFT graphical resource (.jpg, .png, ...) in a figure
% Inputs:
%   rscName:    relative path and name of resource. Path is relative to
%               <sift-root>/resources
%   hax:        optional axis in which to plot image
%
%
% Author: Tim Mullen 2013, SCCN/INC, UCSD.
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


if nargin<2 || isempty(ax)
    hax = [];
end
if isempty(rscName)
    return; end
% remove leading/trailing whitespace
rscName = strtrim(rscName);
if ispc
    strrep(rscName,'/','\');
else
    strrep(rscName,'\','/');
end
% remove any leading slash
if strcmp(rscName(1),filesep)
    rscName(1) = []; end
% get sift rootpath
siftroot = fileparts(utl_whichfile('StartSift.m'));
% construct resource path
rscPath = [siftroot filesep 'resources' filesep rscName];
% check that file exists
if ~exist(rscPath,'file')
    error('SIFT:hlp_viewGraphicsResource:badResourcePath','Resource [%s] does not exist',rscPath);
end
[pathstr name ext] = fileparts(rscPath);
% display resource in figure
im = imread(rscPath);
if isempty(im)
    error('SIFT:hlp_viewGraphicsResource:badResource','Resource [%s] cannot be loaded',rscPath);
end
if isempty(hax)
    figh=figure('Name',[name ext]);
    hax=axes('parent',figh,'position',[0.01 0.01 0.98 0.98]);
end
% render the graphic
image(im,'parent',hax);
axis(hax,'off'); 
axis(hax,'image');