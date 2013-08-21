function cmap = tintmap(tint, M)
%TINTMAP           Tinted monochrome color map.
%   CMAP = TINTMAP(TINT, M), where TINT is a 1-by-3 RGB vector, returns an
%   M-by-3 matrix containing a monochromatic color map based on the color
%   TINT.
%
%   CMAP = TINTMAP(TINT) uses the length of the colormap in the current
%   figure for M.
%   
%   The tinted colormap is made combining the hue and saturation of the
%   color described by TINT with a linear brightness gradient.  A gamma
%   correction of +0.66 is then applied to the resulting rgb colormap.
%
%   See also COLORMAP, GRAY.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 2),  M = size(colormap, 1);  end;
if (length(tint)~=3 || isvectord(tint)~=2 || ~all(tint>=0.0 & tint<=1.0))
	error('TINT must be a 1-by-3 RGB vector with values between 0 and 1.');
end

%%%%%%%%%%%%%%%%%%%%%%%% Add tint via HSV space %%%%%%%%%%%%%%%%%%%%%%
% use hue, sat from TINT w/ gamma-corrected brightness gradient
tint_hsv = repmat(rgb2hsv(tint), M, 1);
tint_hsv(:,3) = linspace(0,1,M)';
cmap = hsv2rgb(tint_hsv).^(2/3);
