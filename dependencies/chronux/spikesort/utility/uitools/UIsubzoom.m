function ax = UIsubzoom(ax)
%UIsubzoom        Adds a context menu for zooming subplots to full figure.
%   UIsubzoom creates a context menu for the current axes that allows the
%   user to move the current axes to a zoomed position filling (most of)
%   the parent figure.  Unchecking the zoom option restores axes to their
%   position before the most recent zoom.
%
%   If the targeted axes contain one or more image objects, the context 
%   menu is added to all image objects in addition to the axes.  If
%   UIsubzoom is called on axes that already have this context menu, it
%   will further link the context menu to any image objects that do not
%   already have it.
%
%   UIsubzoom(HANDLE) associates the menu with the axes specified by
%   HANDLE.
%
%   AXHANDLE = UIsubzoom(...) returns a handle to the parent axes
%   associated with the new context menu.

%state = warning;  warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 1),  ax = gca;  end;
if (isempty(ax)), error('No object currently selected.'); end;
switch get(ax, 'Type'),   % restrict to valid objects ...
	case 'axes',    % do something here
    otherwise,      error(['UIsubzoom not defined for objects of type ' get(h, 'Type') '.']);
end

info.targetaxs = ax;
info.targetfig = get(ax, 'Parent');

%%%%%%%%%%%%%%%%%%%%%%%%% Create Context Menu %%%%%%%%%%%%%%%%%%%%%%%%
imgkids = findobj(ax, 'Type', 'image');  % remember image objs
for h = [imgkids(:)' ax];
	menu = find_uimenu(h, 'subzoom', 'Expand axes', @CB_subzoom);
	if (isempty(get(menu, 'UserData'))),
		set(menu, 'UserData', info);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cleanup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout == 0), clear ax;  end

%warning(state);