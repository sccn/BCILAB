function oldinfo = busyfigure(figurehandle, busy, statusstring, oldinfo)
%BUSYFIGURE        Disable/restore a UI figure.
%   OLDINFO = BUSYFIGURE(FIGUREHANDLE, BUSY) freezes (or unfreezes) the UI
%   figure when it is busy (done) with a task; the BUSY argument is 1 to
%   freeze and 0 to unfreeze.  Freezing disables all uicontrols, disables
%   the figure's click callback, and sets the mouse pointer to a watch.
%
%   BUSYFIGURE(FIGUREHANDLE, BUSY, STATUSSTRING, OLDINFO) also sets status
%   information via the Matlab 'setstatus' function; passing an empty
%   string will skip this call.  (Note that passing in an empty string
%   does _not_ set the status to ''.)  When freezing, the figure's
%   original pre-freeze state is returned; pass this structure as the
%   fourth argument when unfreezing to restore state.  

if (busy)  % FREEZE

	% Get current uicontrol handles and enable states.
    controls = get(figurehandle, 'Children');
    oldinfo.controls = controls(strcmp(get(controls, 'Type'), 'uicontrol'));
	oldinfo.enables = get(oldinfo.controls, 'Enable');
	set(oldinfo.controls, 'Enable', 'off');

	% Deal with the status bar.
	statusbarH = controls(find(strcmp(get(controls, 'Tag'), 'Status')));
    if (~isempty(statusbarH)), set(statusbarH, 'Enable', 'on'); end   % need enabled
    if (~isempty(statusstring)), setstatus(figurehandle, statusstring); end
	
	% Unset the figure mouse press callback, but save it for later restoration.
    oldinfo.figure_click_callback = get(figurehandle, 'ButtonDownFcn');
    set(figurehandle, 'ButtonDownFcn', '');
    
	% Set pointer to watch to visually indicate that we're busy.
    oldinfo.pointer = get(figurehandle, 'Pointer');
    setptr(figurehandle, 'watch');
	
else       % RESTORE
	if (nargin < 4)
		error('Unfreezing requires the structure returned when the figure was frozen.');
	end

	% Restore enable states.
    for cntrl = 1:length(oldinfo.enables)
        set(oldinfo.controls(cntrl), 'Enable', oldinfo.enables{cntrl});
    end

	% Status update.
	if (~isempty(statusstring)), setstatus(figurehandle, statusstring); end

	% Reset the figure mouse click callback.
    set(figurehandle, 'ButtonDownFcn', oldinfo.figure_click_callback);
    
	% Finally, restore old pointer.
	setptr(figurehandle, oldinfo.pointer);
end
