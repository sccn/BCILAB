function [elbowVal elbowIdx] = hlp_findElbow(curve)
% [elbowVal elbowIdx] = hlp_findElbow(curve)
%
% Find the 'elbow' in a curve. We do this by drawing a line from the first 
% to the last point of the curve and then finding the data point that is 
% farthest away from that line. Technically, if b is the vector from the 
% origin (first value in curve) to the last value in curve, for each point 
% p on the curve, we find the one with the maximum distance d given by:
%
% d     = |p-(p.b_hat)b_hat|
% b_hat = b/|b|
%
% adapted from [2]
%
% See Also: est_selModelOrder()
%
% References:
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT
% [2] http://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
%
% Author: @Amro 01/07/2011, Tim Mullen 04/26/2011, (c) SCCN/INC, UCSD. 
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

% make row vector
curve = curve(:)';

% get coordinates of all the points
nPoints = length(curve);
allCoord = [1:nPoints;curve]';

% pull out first point
firstPoint = allCoord(1,:);

% get vector between first and last point - this is the line
lineVec = allCoord(end,:) - firstPoint;

% normalize the line vector
lineVecN = lineVec / sqrt(sum(lineVec.^2));

% find the distance from each point to the line:
% vector between all points and first point
vecFromFirst = bsxfun(@minus, allCoord, firstPoint);

% To calculate the distance to the line, we split vecFromFirst into two 
% components, one that is parallel to the line and one that is perpendicular 
% Then, we take the norm of the part that is perpendicular to the line and 
% get the distance.
% We find the vector parallel to the line by projecting vecFromFirst onto 
% the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
% We project vecFromFirst by taking the scalar product of the vector with 
% the unit vector that points in the direction of the line (this gives us 
% the length of the projection of vecFromFirst onto the line). If we 
% multiply the scalar product by the unit vector, we have vecFromFirstParallel
scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1), 2);
vecFromFirstParallel = scalarProduct * lineVecN;
vecToLine = vecFromFirst - vecFromFirstParallel;

% distance to line is the norm of vecToLine
distToLine = sqrt(sum(vecToLine.^2,2));

% now all you need is to find the maximum
[dummy,elbowIdx] = max(distToLine);
elbowVal = curve(elbowIdx);

