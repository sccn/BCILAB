function h = UIcolorshift(varargin)
%UIcolorshift      Creates an interactive colorbar.
%   UIcolorshift adds a (vertical) colorbar to the current axes.  This
%   colorbar can be used to shift or rescale the color mapping in the
%   current axes.  Clicking near the top of the colorbar axes
%   up or down (right or left for horizontal colorbars) changes the value
%   mapped to the top color scale limit. Similarly, clicking near the
%   bottom changes the value mapped to the bottom color scale limit.
%
%   Clicking elsewhere and dragging cyclically shifts the color map.  This
%   shift is applied to the whole figure.
%
%   UIcolorshift(ARG1,...) takes the same arguments as COLORBAR.
%
%   H = UIcolorshift(...) returns a handle to the colorbar axes.

%%%%%%%%%%%%%%%%%%%%%%%%% Create the Colorbar %%%%%%%%%%%%%%%%%%%%%%%%
target = gca;
fighdl = gcf;
h = colorbar(varargin{:});  % create or refresh an existing colorbar
%h = find_colorbar(target);  % <=R13 only

%%%%%%%%%%%%%%%%%%%%% Get a handle to the Colorbar %%%%%%%%%%%%%%%%%%%
cimage = findobj(h, 'Type', 'image');

cbar_pos = get(h, 'Position');
ax_pos   = get(target, 'Position');
if (cbar_pos(1) > sum(ax_pos([1,3]))),   mode = 'V';
else                                     mode = 'H';
end

%%%%%%%%%%%%%%%%%%%%%%%%% Assign the Callback %%%%%%%%%%%%%%%%%%%%%%%%
userdata = get(h,'UserData');
userdata.fighdl = fighdl;
userdata.cbrhdl = h;
userdata.target = target;
userdata.cimage = cimage;
userdata.mode   = mode;
set(h, 'UserData', userdata, 'ButtonDownFcn', @CB_colorshift);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cleanup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout == 0), clear h;  end
axes(target);   % reset current axes