function CB_movieplayer(handle, events)
%CB_MOVIEPLAYER    Callback for UImyfunc.

persistent running direction renderer;
if (isempty(running)), running = 0;  end;
if (isempty(direction)), direction = 1;  end;

info  = guidata(handle);  % get data & associated info

% Constants
fpsupdate = 8;

switch(handle),
	case {info.menu.intp, info.menu.rndr}       % TOGGLE COLOR INTERPOLATION
		% OpenGL is sometimes faster, depending on hardware.
		% However, OpenGL interpolates in RGB rather than through the colormap.
		togglecheck(handle);
		if (onoff2bool(get(info.menu.intp,'Checked')) || info.hires),
			shade = 'interp';  xlim = [1,info.N];   ylim = [1,info.M];   % don't need pad
		else
	        shade = 'flat';    xlim = [1,info.N+1]; ylim = [1,info.M+1]; % need pad'd row/col
		end
		if (onoff2bool(get(info.menu.rndr, 'Checked'))),
			renderer = 'opengl';
		else
			renderer = 'zbuffer';
		end
		shading(info.haxs, shade);   set(info.hfig, 'Renderer', renderer);
		set(info.haxs, 'XLim', xlim, 'YLim', ylim);	
		
	case ({info.menu.play, info.menu.rply}),   % FORWARD OR REVERSE PLAY
		check = togglecheck(handle);
		renderer = get(info.hfig, 'Renderer');		
		if (check)  % starting to play
			if (handle == info.menu.play)  % prevent switching directions while playing
				direction =  1;  set(info.menu.rply, 'Enable', 'off');
            else
				direction = -1;  set(info.menu.play, 'Enable', 'off');
			end;
			set(info.menu.parm, 'Enable', 'off');
			running = 1;   % anyone else can set this to 0 to interrupt play
			counter = 1;   fps = 0;   tic;   % keep track of display frame rate
			next = info.frame;
			while (running)
				next = next + direction*(info.skip);   % try to increment ...
				if (next>info.P), direction = -1;  break;  end;
				if (next < 1),    direction =  1;  break; end;
				updateframe(info, next);   info.frame = next;
				set(info.hfig, 'Name', sprintf('%s  FPS: % 5.3f', renderer, fps));
				drawnow;  pause(info.slow);

				if (counter == fpsupdate)   % frame rate updates
					counter = 1;   fps = fpsupdate/toc;  tic;
                else  counter = counter + 1;
				end
			end
			running = 0;
			set([info.menu.play, info.menu.rply, info.menu.parm], 'Enable', 'on', 'Checked', 'off');
			set(info.hfig, 'Name', 'UImovieplayer');
			guidata(info.hfig, info);
		else
			running = 0;   % signal an interruption
		end
		
	case info.menu.parm,     % TIMING PARAMETERS
		params = {'SKIP: frames advanced per step', ...
				  'SLOW: pause interval during playback (msec)', ...
				  'Fs:   temporal sampling rate (Hz)'};
		defaults = {sprintf('%d', info.skip), sprintf('%5.3f', info.slow), sprintf('%8.6f', info.Fs)};
		values = inputdlg(params, '', 1, defaults);
		if (~isempty(values))   % validate parameters
			skip = str2num(values{1});  skip = round(skip);
			if (~isempty(skip) && skip > 0 && skip < info.P), info.skip = skip;  end;
			slow = str2num(values{2});
			if (~isempty(slow) && slow > 0), info.slow = slow;  end;
			Fs   = str2num(values{3});
			if (~isempty(Fs) && Fs > 0), info.T = (1:info.P)/Fs;  info.Fs = Fs;  end;
		end
		updateframe(info, info.frame);
		info = redraw_timeseries(info);
		guidata(info.hfig, info);
		
	case info.hfig,          % KEY PRESS HANDLER
		key = get(info.hfig, 'CurrentChar');
		if (key == ' ')   % pause always works
			if (running),  running = 0;
            else 
				if (direction == 1),  CB_movieplayer(info.menu.play);
                else                  CB_movieplayer(info.menu.rply);
				end
			end
		end
		if (~running)     % but other keypresses are ignored if movie is playing 
			next = info.frame;
			switch(key)
				case '<',  next = 1;
				case ',',  next = next - info.skip;
				case '?',  % next = next;
				case '.',  next = next + info.skip;
				case '>',  next = info.P;
				otherwise, return;
			end
			if (next > 0 && next <= info.P)
				updateframe(info, next);   info.frame = next;
			end
			guidata(info.hfig, info);
		end
		
	case info.hdat,         % TIME SERIES HANDLER
		% clicking in the image draws a time series for corresponding
		% pixel in the lower time-series axes
		if (~running && strcmp(get(info.hfig, 'SelectionType'), 'normal'))
 			click = floor(get(info.haxs, 'CurrentPoint'));  click = click(1,1:2)';
			if (click(1)>=1 && click(1)<=info.N && click(2)>=1 && click(2)<=info.M)
				info.tslist = unique([click info.tslist]', 'rows')';
				info = redraw_timeseries(info);
				guidata(info.hfig, info);
			end
		end
		
	case info.hax2,         % TIME SERIES AXES
		if (~running)
			switch(get(info.hfig, 'SelectionType')),
				case 'open',  % dbl-click clears the axes
					cla;   legend off;
					if (ishandle(info.hdtx)),  delete(info.hdtx);  end;
					info.tslist = [];   info.hcsr = [];   info.hdtx = [];
					guidata(info.hfig, info);
				otherwise,
			end
		end
	
	case {info.menu.dxdt, info.menu.legd}    % TIME SERIES AXES -- special functions
		if (~running)
			check = togglecheck(handle);
			info = redraw_timeseries(info);
			guidata(info.hfig, info);
		end
		
	case info.hcsr,
		if (~running && strcmp(get(info.hfig,'SelectionType'),'normal'))
			wide = get(info.hax2, 'Position');  wide = wide(3);				
			xlim = get(info.hax2, 'XLim');
			dataperpixel = (xlim(2)-xlim(1)) ./ wide;  % cursor motion scale
			
			time = get(handle, 'XData');    time = time(1);  % starting pt
			
			info.buttondown = 1;   guidata(handle, info);
			pntA = get(info.hfig, 'CurrentPoint');
			while (info.buttondown)
				running = 1;  % make sure that we're not interrupted
				pntB = get(info.hfig, 'CurrentPoint');
				newtime = time + (pntB(1)-pntA(1))*dataperpixel;  
				newtime = max(min(newtime, xlim(2)), xlim(1));  % clip to axes
				newframe = round(newtime*info.Fs);
				newframe = max(min(newframe,info.P),1);

				set(info.hcsr, 'XData',  [newtime newtime]);
				updateframe(info, newframe);     drawnow;
				info = guidata(handle);
			end
			running = 0;  info.frame = newframe;
			guidata(handle, info);
		end
		
	otherwise,
		error('Error in UImovieplayer callback: unknown function.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% HELPER Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check = togglecheck(handle)
% When HANDLE is a handle to a UI menu, its 'Checked' property
%  is toggled between 'off' and 'on', and a boolean value representing
%  the new state ('off'->0, 'on'->1) is returned.
check = ~onoff2bool(get(handle,'Checked'));
set(handle, 'Checked', bool2onoff(check));

function updateframe(info, framenum)
% Updates data and title/time displays
set(info.httl, 'String', sprintf('Time: % 8.3f', info.T(framenum)));
set(info.hdat, 'CData',  [double(info.data(:,:,framenum)) info.COLpad; info.ROWpad]);
if (~isempty(info.hcsr))
	set(info.hcsr, 'XData',  [info.T(framenum) info.T(framenum)]);
end


function info = redraw_timeseries(info)
% Redraw the time-series axes
if (~isempty(info.tslist))
	N = size(info.tslist,2);
	
	axes(info.hax2);
	subdata = [];  leg = {};
	for k = 1:N    % extract pixel time series; do this fresh each time
		subdata = cat(2, subdata, squeeze(info.data(info.tslist(2,k), info.tslist(1,k), :)));
		leg = cat(2, leg, {sprintf('R:%d/C:%d',info.tslist(2,k),info.tslist(1,k))}); % legend
	end
	
	clim = info.CLIM;
	if (onoff2bool(get(info.menu.dxdt, 'Checked'))),
        subdata = diff(double([subdata(1,:); subdata]), [], 1);
		clim = minmax(subdata);
	end
	axes(info.hax2);  plot(info.T, subdata);  set(info.hax2, 'XLim', [info.T(1), info.T(end)]);
	if (onoff2bool(get(info.menu.legd, 'Checked'))), legend(leg);  end;
	
	% reset the button down callback & draw reference line
	info.hcsr = line([0 0], clim, 'LineWidth', 2, 'Color', [0.6 0.6 0.6], 'LineStyle', '--');
	set([info.hcsr,info.hax2], 'ButtonDownFcn', @CB_movieplayer);
	
	% Draw markers on the data to show time series source points
	axes(info.haxs);  hold on;
	if (ishandle(info.hdtx)),  delete(info.hdtx);  end;
	pts = cat(1, mat2cell(info.tslist(1,:),1, ones(N,1)), ...  % dirty trick to plot pts in diff colors
		       mat2cell(info.tslist(2,:),1,ones(N,1)), mat2cell(repmat('x',1,N),1,ones(N,1)));
	hold on;  info.hdtx = plot(pts{:});  hold off;
	set(info.hax2, 'YLim', clim);
end
