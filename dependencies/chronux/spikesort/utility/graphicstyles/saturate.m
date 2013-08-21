function newmap = saturate(cmap)
%SATURATE          Add saturation markers to a colormap.
%   SATURATE with no arguments modifies the color map of the current
%   figure to emphasize minimal and maximal values.  The actual change
%   is chosen heuristically to contrast with the current color map.  By
%   default, the maximum value is set to white and the minumum value is
%   set to black.  If either extremal value in the existing colormap is
%   already close (i.e., within a radius of 0.02 in RGB-space) to black or
%   white, red and blue are used for the max and min values, respectively.
%
%   NEWMAP = SATURATE retrieves the colormap from the current figure,
%   adds saturation markers and returns the result without modifying the
%   current figure. 
%
%   NEWMAP = SATURATE(MAP) adds saturation markers to the colormap MAP and
%   returns the altered colormap.

bwcutoff = 0.02;

%%%%%  Argument checking.
if (nargin < 1)
	cmap = colormap;
elseif ((size(cmap,2) ~= 3) || ~all(cmap(:)>=0) || ~all(cmap(:)<=1))
	error('Colormap must have 3 columns: [R,G,B] with all values in [0,1].');
end

%%%%%%%%%%%%%%%%%% Decide what modification to make %%%%%%%%%%%%%%%%%%
over = [1 1 1];  under = [0 0 0];

if    (norm(cmap(1,:)-[0 0 0])<bwcutoff   || norm(cmap(1,:)-[1 1 1])<bwcutoff ...
    || norm(cmap(end,:)-[0 0 0])<bwcutoff || norm(cmap(end,:)-[1 1 1])<bwcutoff)
	over = [1 0 0];  under = [0 0 1];
end

%%%%%%%%%%%%%%%%%%%% Set the colormap or return it %%%%%%%%%%%%%%%%%%%%
cmap(1,:) = under;   cmap(end,:) = over;
if (nargout == 0),  colormap(cmap);
else                newmap = cmap;
end

