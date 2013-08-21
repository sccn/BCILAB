function CB_colorshift(handle, events)
%CB_COLORSHIFT     Callback for UIcolorshift.

% Fields for UserData structure:
%  fighdl - parent figure hndl    mouseroot - global coords initial click
%  target - data axes hndl        climsinit - initial clim
%  cimage - colorbar image hndl   limitSEL  - 1 = bottom limit, 2 = top
%  mode   - 'H'oriz or 'V'ert     mouseSEL  - 1 = horiz (X), 2 = vert (Y)
%  cbrhdl - colorbar axes hndl    limitstr  - 'XLim'=horiz, 'YLim'=vert
%  cmap   - initial colormap      cmaplen - # rows in colormap
%  datastr - 'XData'=horiz,'YData'=vert
%  incr_per_pix - clim incr per pixel of mouse motion

%  shft_per_pix - clim shift per pixel of mouse motion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set Up Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
info = get(handle, 'UserData');

info.mouseSEL = (info.mode == 'V') + 1;   % 1 for 'H' mode, 2 for 'V'
if (info.mode == 'H'),  info.limitstr = 'XLim';  info.datastr = 'XData';
else                    info.limitstr = 'YLim';  info.datastr = 'YData';
end

info.climsinit = get(info.target,'CLim');   
info.incr_per_pix = diff(info.climsinit)./200;

info.cmap = get(info.fighdl, 'Colormap');
info.cmaplen = size(info.cmap,1);
info.shft_per_pix = info.cmaplen./256;


% mouseclick in data units
clickData = get(handle,'CurrentPoint');  clickData = clickData(1,info.mouseSEL);

% mouseclick in figure units
clickRoot = get(0, 'PointerLocation');   info.mouseinit = clickRoot(info.mouseSEL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Decide Action %%%%%%%%%%%%%%%%%%%%%%%%%%%
dlims = get(handle, info.limitstr);              
click = (clickData - dlims(1)) ./ (dlims(2)-dlims(1));  % data units -> normalized
if (click >= 0.95),     mode = 'top';   info.limitSEL = 2;
elseif (click <= 0.05), mode = 'bot';   info.limitSEL = 1;
else                    mode = 'shf';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Track Mouse %%%%%%%%%%%%%%%%%%%%%%%%%%%%
oldMotionFcn   = get(info.fighdl, 'WindowButtonMotionFcn');
oldButtonUpFcn = get(info.fighdl, 'WindowButtonUpFcn');
set(info.fighdl, 'WindowButtonUpFcn', {@done_cbar});
switch(mode),
    case {'top', 'bot'}
        set(info.fighdl, 'WindowButtonMotionFcn', {@scale_cbar, info});
    case 'shf',
        set(info.fighdl, 'WindowButtonMotionFcn', {@shift_cbar, info});
end
waitfor(info.fighdl, 'WindowButtonMotionFcn', '');

set(info.fighdl, 'WindowButtonMotionFcn', oldMotionFcn);
set(info.fighdl, 'WindowButtonUpFcn', oldButtonUpFcn);

axes(info.target);  % R13 TMW bug on deletion if colorbar remains current axes
if (info.mode == 'H')
    UIcolorshift('horiz', 'peer',info.target);   % refresh display
else
    UIcolorshift('peer', info.target);  
end
% click = round((click-range(1)) ./ (range(2)-range(1)) .* 255);
% shift = mod([-127:128] + click, 256) + 1;      % make shifted indices
% cmap = colormap; 
% colormap(cmap(shift,:));  % then shift the map


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Callbacks for the Callback %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Used to signal to cb_colorshift that tracking is over.
function done_cbar(handle, event)
set(handle, 'WindowButtonMotionFcn', '')

%%%%% Rescale color limits
function scale_cbar(handle, event, ud)
mouse = get(0,'PointerLocation');
ud.climsinit(ud.limitSEL) = ud.climsinit(ud.limitSEL) - (mouse(ud.mouseSEL)-ud.mouseinit)*ud.incr_per_pix;
if (diff(ud.climsinit) > 0)
    set(ud.target, 'CLim', ud.climsinit);
    set(ud.cbrhdl, ud.limitstr, ud.climsinit);
    set(ud.cimage, ud.datastr, ud.climsinit);
end

%%%%% Shift color limits
function shift_cbar(handle, event, ud)
mouse = get(0,'PointerLocation');
shift = round((mouse(ud.mouseSEL)-ud.mouseinit)*ud.shft_per_pix);
cmap = ud.cmap;
set(ud.fighdl, 'Colormap', cmap(circshift((1:ud.cmaplen)', shift), :));