function hfig = UImovieplayer(data, Fs, frame)
%UImovieplayer        Show a 3-D matrix as an XY-T movie.
%   UImovieplayer(DATA), for a M x N x P matrix DATA, creates a new figure
%   and shows the frames DATA(:,:,k), 1 <= k <= P, as a movie.
% 
%   UImovieplayer(DATA, Fs) specifies the temporal sampling rate of the
%   frames.  This value is used for displaying time and for temporal
%   filtering on the time-series axes (see below).
%
%   UImovieplayer(DATA, Fs, FRAME), where 1 <= FRAME <= P, specifies the
%   index of the frame to display when the UImovieplayer starts.
%  
%   Several tools control the resulting movieplayer figure:
%   
%   The following keystroke commands change frames (all except Pause are
%   disabled during automatic playback):
%       '<'  Rewind to first frame (frame 1)
%       ','  Rewind by 1 frame   (or 'Skip' frames, see above)
%       ' '  (space bar)  Pause/unpause automatic playback.
%       '.'  Advance by 1 frame  (or 'Skip' frames, see above)
%       '>'  Advance to last frame (frame P)
%
%   Right clicking the data invokes a menu with the following commands: 
%      Play          Automatic playback (frame 1 -> frame P).
%      Reverse       Automatic playback in reverse (frame P -> frame 1).
%      Interpolate   Toggles interpolated color display.
%      OpenGL        Toggles between the OpenGL and zbuffer renderers.
%      Params        Creates a dialog box to ask for timing parameters:
%         Skip  (default: 1 frms) # frames to advance btw frames.
%         Slow  (default: 0 msec) # msec to pause btw frame advances.
%         Fs    (default: 1 Hz  ) temporal sampling rate
%
%   Left clicking over the data plots the time-series for the clicked
%   pixel in the lower axes.  Double click on the time-series axes
%   background to clear them.
%
%   Left click and drag the time cursor on the time-series axes to rapidly
%   move to another part of the recording.
%
%   Right-clicking the time-series axes invokes these menu options:
%     Expand Axes       Toggles a zoomed view; see UIsubzoom.
%     Temporal Diff     Temporal differences of time-series data.
%     Show Legend       Toggles legend display.
%
%   The colorbar alters the color scale; see UIcolorshift.
%
%   H = UImovieplayer(...) returns a handle to the new movieplayer figure.
%
%   Technical note:
%   The interpolation option is disabled when the number of spatial data
%   points (M x N) exceeds 2500 (since such large movies are rendered
%   using a technique that does not support interpolation).  However, when
%   M x N is less than 2500, interpolation can be used for spatial data
%   smoothing.  Note that OpenGL and zbuffer interpolate differently; in
%   general, zbuffer interpolates in the colormap and is thus a more
%   accurate spatial smoothing.  However, OpenGL may be faster on systems
%   where hardware OpenGL acceleration is available.

% programmer note -- at some point, should allow non-square X-Y grid 

%state = warning;  warning off;
surfimg_cutoff = 2500;    % pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if     (nargin == 1),  frame = 1;   Fs = 1;
elseif (nargin == 2),  frame = 1;
end

[M,N,P] = size(data);    T = (1:P)./Fs;    AVG = mean(data,3);
if (M*N < surfimg_cutoff),  hires = 0;  else  hires = 1;  end;
ROWpad = NaN*ones(1,N+1);    COLpad = NaN*ones(M,1);

slow = 0;    skip = 1;   tslist = [];
if (frame<1 || frame>P),  error('Requested start frame does not exist.');  end;
if (Fs<=0),               error('Sampling rate must be positive.'); end;

if (~isreal(data)),  error('DATA must be real');   % get data limits
elseif  (isa(data, 'uint8')), CLIM =  [   0 255];
elseif  (isa(data, 'int8')),  CLIM =  [-128 127];
elseif  (isa(data, 'double')), CLIM = minmax(data);  %    (faster then min/max)
else  CLIM(1) = min(data(:));  CLIM(2) = max(data(:));
end

hfig = figure('Position', [500 300 600 600]);
subplot(5,1,5);    hax2 = gca;    hcsr = [];  % cursor invisible
subplot(5,1,1:4);  haxs = gca;    hdtx = [];  % no time series to start ...
if (hires),  hdat = imagesc(data(:,:,frame));
else         hdat = pcolor(1:(N+1),1:(M+1),[double(data(:,:,1)) COLpad; ROWpad]);
end
saturate;
httl = title('');
hxlb =  ylabel('See UImovieplayer help for controls.');

%%%%%%%%%%%%%%%%%%% Movieplayer Control Callbacks %%%%%%%%%%%%%%%%%%%%%
cmenu = uicontextmenu;   set(hdat, 'UIContextMenu', cmenu);

menu.play = uimenu(cmenu, 'Checked', 'off', 'Label', 'Play', 'Callback', ...
	              @CB_movieplayer);
menu.rply = uimenu(cmenu, 'Checked', 'off', 'Label', 'Reverse', 'Callback', ...
	              @CB_movieplayer);

menu.intp = uimenu(cmenu, 'Checked', 'off',  'Label', 'Interpolate', 'Callback', ...,
	              @CB_movieplayer, 'Separator', 'on');
menu.rndr = uimenu(cmenu, 'Checked', 'on',  'Label', 'OpenGL', 'Callback', ...,
	              @CB_movieplayer);  % toggled on->off by default below
 
menu.parm = uimenu(cmenu, 'Checked', 'off', 'Label', 'Timing Params', 'Callback', ...,
	              @CB_movieplayer, 'Separator', 'on');

if (hires), set(menu.intp, 'Checked', 'off', 'Enable', 'off'); end;
set(hfig, 'KeyPressFcn', @CB_movieplayer);
set(hdat, 'ButtonDownFcn', @CB_movieplayer);

set(hax2, 'ButtonDownFcn', @CB_movieplayer, 'Units', 'pixels');  UIsubzoom(hax2);
cmenu2 = get(hax2, 'UIContextMenu');
menu.dxdt = uimenu(cmenu2, 'Tag', 'dx/dt', 'Checked', 'off', ...
			       'Label', 'Temporal Diff', 'Callback', @CB_movieplayer);
menu.legd = uimenu(cmenu2, 'Tag', 'ShowLegend', 'Checked', 'off', ...
			       'Label', 'Show Legend', 'Callback', @CB_movieplayer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Store UserData %%%%%%%%%%%%%%%%%%%%%%%%%%
udlist = {'data', 'N', 'M', 'P', 'T', 'CLIM', 'AVG', 'hires',  ...
		  'hdat', 'hfig', 'haxs', 'hax2', 'httl', 'hcsr', 'hdtx', 'cmenu', 'cmenu2', 'menu', ...
	  	  'frame', 'slow', 'skip', 'Fs', 'tslist', 'COLpad', 'ROWpad',};
guidata(hfig, structpack(udlist));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Default State %%%%%%%%%%%%%%%%%%%%%%%%%%
% Signaling function for mouse presses
set(hfig, 'WindowButtonUpFcn', 'h = guidata(gcbo);  h.buttondown = 0;  guidata(gcbo, h);');
set(hfig, 'WindowButtonMotionFcn', '% dummy function to force CurrentPoint updates');
set(hfig, 'Name', 'UImovieplayer', 'NumberTitle', 'off');
set(haxs, 'CLim', CLIM);  UIcolorshift;
axis xy square equal;  % default since we're usu looking at data instead of an image
CB_movieplayer(menu.rndr);  % force a redraw to set shading/limits
set(hfig, 'CurrentCharacter', '?');  CB_movieplayer(hfig);  % force redraw for text

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cleanup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout == 0), clear hfig;  end
%warning(state);
