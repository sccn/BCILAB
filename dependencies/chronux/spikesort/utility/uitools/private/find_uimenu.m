function [cxitem,cxmenu] = find_uimenu(parent, tag, label, callback)
%FIND_UIMENU        Finds/creates an item in a UIcontextmenu.
%   HITEM = FIND_UIMENU(PARENT,TAG) searches the UIcontextmenu (if any)
%   associated with the object PARENT for a menu item with tag TAG.  If
%   such a menu item exists, its handle is returned in HITEM.  If no such
%   item is found, or if PARENT does not have an associated UIcontextmenu,
%   the function returns the empty matrix.
%
%   HITEM = FIND_UIMENU(PARENT,TAG,LABEL,CALLBACK) similarly searches the
%   PARENT UIcontextmenu for an item with tag TAG.  If no such menu is
%   found, it is created (along with a UIcontextmenu, if necessary) with
%   the following properties:  'Tag' is set to TAG, 'Label' is set to
%   LABEL, and 'Checked' is set to 'off'.
%
%   [HITEM,HMENU] = FIND_UIMENU(...) also returns a handle to the
%   UIcontextmenu itself.


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~ishandle(parent)),  error('First argument must be a valid handle.');  end;
searchonly = (nargin < 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Find/create Menu %%%%%%%%%%%%%%%%%%%%%%%%%
cxmenu = get(parent, 'UIContextMenu');
if (isempty(cxmenu)),   % if no context menu exists ...
	if (searchonly),  cxitem = [];   return;
    else              cxmenu = uicontextmenu; set(parent, 'UIContextMenu', cxmenu);
	end
end

cxitem = findobj(cxmenu, 'Tag', tag);
if (isempty(cxitem))    % if no matching menu item exists ...
	if (searchonly),  return;
    else              cxitem = uimenu(cxmenu, 'Tag', tag, 'Checked', 'off', ...
			                          'Label', label, 'Callback', callback);
	end
end
