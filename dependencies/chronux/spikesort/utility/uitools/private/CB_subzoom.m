function CB_subzoom(handle, event)
%CB_MYFUNC         Callback for UImyfunc.

hei_zoom = 0.85;  % fractional height in zoomed state
wid_zoom = 0.85;  % fractional width in zoomed state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Access Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First we get the target axes & its context menu -- these are in charge.
userdata = get(handle, 'UserData');  
target = userdata.targetaxs;
mymenu = find_uimenu(target, 'subzoom');
% If any kids have similar menu's, they need to mirror state
allcontexts = get(get(target, 'Children'), 'UIContextMenu');
if (iscell(allcontexts)),  allcontexts = cat(1, allcontexts{:});  end;
chmenu = findobj(allcontexts, 'Tag', 'subzoom');

userdata = get(mymenu, 'UserData');

userdata.zoomstate = ~onoff2bool(get(mymenu,'Checked'));
set([mymenu, chmenu], 'Checked', bool2onoff(userdata.zoomstate));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oldunits = get(target, 'Units');   set(target, 'Units', 'normalized');
if (userdata.zoomstate)   % ZOOMING IN
	userdata.revert = get(target, 'Position');
	set(target, 'Position', ...
	             [(1-wid_zoom)/2 (1-hei_zoom)/2 wid_zoom hei_zoom]);
	uistack(target, 'top');
else                      % ZOOMING OUT
	set(target, 'Position', userdata.revert);
	userdata.revert = [];
end
if (~isempty(legend(target))),  legend(target);  end;
set(target, 'Units', oldunits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clean Up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(mymenu, 'UserData', userdata);