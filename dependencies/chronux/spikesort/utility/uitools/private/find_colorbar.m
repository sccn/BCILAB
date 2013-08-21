function cbarax = find_colorbar(ax)
%FIND_COLORBAR     Finds the colorbar associated with a given axes.
%   CBARAX = FIND_COLORBAR(AX) returns the handle of the colorbar
%   associated with the axes AX.  If no colorbar is found, CBAR_AX will be
%   empty.

%%%%% TMW's COLORBAR function creates a hidden proxy text object in a 
%%%%% plot when it associates a colorbar with it.
shh = get(0, 'ShowHiddenHandles');   set(0, 'ShowHiddenHandles', 'on');
proxy = findobj(get(ax, 'Children'), 'Tag', 'ColorbarDeleteProxy');
set(0, 'ShowHiddenHandles', shh);

%%%%% The proxy's UserData contains a handle to the colorbar ...
cbarax = get(proxy, 'UserData');  % (get will return [] if proxy is [])