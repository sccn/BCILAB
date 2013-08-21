function outdata=specscope(indata)
% record and plot audio spectrogram
%
% Usage: outdata=specscope(indata)
%
%  Input: indata (optional)
%    Displays a recorded piece of data, if an argument is passed
%    Otherwise displays audio data from an attached microphone
%
%  Output: outdata (optional)
%    If present, will return up to 10 minutes
%    of captured audio data.
%
%  Note: Parameters such as sampling frequency, number of tapers
%    and display refresh rate may be set below if desired.
%    You can also acquire data from a national instruments card.
%

close all;

%=======================================================================
% Check toolboxes
%=======================================================================

% Check for toolboxes
if not(exist('analoginput','file'));
    fprintf('You need to install the DAQ toolbox first\n');
    return
end
if not(exist('mtspecgramc','file'));
    fprintf('You need to install the Chronux toolbox first from http://chronux.org/\n');
    return
end

%=======================================================================
% Set parameters
%=======================================================================

global acq;

% Set defaults
acq.params.Fs = 44100;
acq.pause = 0;
acq.skips = 0;
acq.stop = 0;
acq.restart = 0;
acq.plot_frequency = 10;
acq.samples_acquired = 0;
acq.spectra = [];
acq.times = [];
defaults
audio_instr;
fig=create_ui;

%=======================================================================
% Check arguments, start DAQ
%=======================================================================

if nargout 
   % save up to ten minutes data, preallocated...
    fprintf('Pre-allocating memory for data save - please be patient.\n');
   outdata=zeros( (acq.params.Fs * 60 * 10), 1 ); 
end    

if nargin == 1;
    acq.indata = indata;
    acq.live = 0;
else
    % Create and set up and start analog input daq
    input=1;
    if input==1;
        acq.ai = analoginput('winsound');
        addchannel( acq.ai, 1 );
    else
        acq.ai = analoginput('nidaq', 1);
        addchannel(acq.ai, 0);
        set(acq.ai,'InputType','SingleEnded')
        set(acq.ai,'TransferMode','Interrupts')
        set(acq.ai,'TriggerType','Manual');
    end
    set( acq.ai, 'SampleRate', acq.params.Fs )
    acq.params.Fs = get( acq.ai, 'SampleRate' );
    set( acq.ai, 'SamplesPerTrigger', inf )
    start(acq.ai)
    acq.live = 1;

    if input==2;
        trigger(acq.ai);
    end
end

acq.samples_per_frame = acq.params.Fs / acq.plot_frequency;
 

%=======================================================================
% The scope main loop
%=======================================================================

acq.t0=clock;
acq.tn=clock;

% Loop over frames to acquire and display
while 1;

    % Check for quit signal
    if acq.stop;
        break;
    end
    
    % Calculate times
    calctime = acq.samples_acquired / acq.params.Fs;
    acq.samples_acquired = acq.samples_acquired + acq.samples_per_frame;
    acq.t1 = clock;
    elapsed = etime(acq.t1,acq.t0);

    % Get a small snippet of data
    if ( acq.live )
       data = getdata( acq.ai, acq.samples_per_frame );
    else
      while elapsed < acq.samples_acquired / acq.params.Fs
        pause( acq.samples_acquired / (acq.params.Fs) - elapsed );
        acq.t1=clock;
        elapsed = etime(acq.t1,acq.t0);
      end
      if acq.samples_acquired + 2 * acq.samples_per_frame >= length( acq.indata )
        acq.stop=1;
      end
      data = acq.indata(floor(acq.samples_acquired+1):floor(acq.samples_per_frame+acq.samples_acquired));
    end

    if nargout 
        outdata(floor(acq.samples_acquired+1):floor(acq.samples_acquired+length(data))) = data(:);
    end

    if acq.restart;
       acq.restart = 0;
       acq.spectra = [];
       acq.times = [];
    end

    % Calculate spectrogram of data snippet
    if acq.deriv
      [s, t, f] = mtspecgramc(diff(data), acq.moving_window, acq.params );
    else
      [s, t, f] = mtspecgramc(data, acq.moving_window, acq.params );
    end
    
    % Add new spectra to that already calculated
    acq.times = [acq.times t+calctime];
    if acq.log
        acq.spectra = [acq.spectra log(s')];
    else
        acq.spectra = [acq.spectra s'];
    end			    

    % Remove old spectra once window reaches desired size
    while acq.times(1,end) - acq.times(1,1) > acq.display_size;

        % Ring buffer!
        y = length(t);
        acq.times(:,1:y) = [];
        acq.spectra(:,1:y) = [];
        
    end

    % Only plot if display is keeping up with real time and not paused
    show_plot=1;
    if nargin==0
       if get(acq.ai, 'SamplesAvailable' ) > 10 * acq.samples_per_frame && acq.pause==0
	  show_plot=0;
       end
    else
      if elapsed > calctime + 0.5
	show_plot=0;
      end
    end

    if acq.pause
       show_plot=0;
    end
    if show_plot
        
        if  acq.bgsub
            acq.mean_spectra = mean( acq.spectra, 2 );
        end
        
        % Normalize until full screen passes by if requested
        if acq.normalize>=1;
            if acq.normalize==1
                acq.tn=clock;
                acq.normalize=2;
            end
            if etime(clock,acq.tn)>1.25*acq.display_size
                acq.normalize=0;
            end
            mins = min(min(acq.spectra));
            maxs = max(max(acq.spectra));
        end

        % Scale the spectra based upon current offset and scale
        if acq.bgsub
           scaled_spectra = acq.offset + ( acq.scale ) * ( acq.spectra - repmat( acq.mean_spectra, [1,size(acq.spectra,2)]) ) / ( maxs - mins ); 
        else
            scaled_spectra = acq.offset + acq.scale * ( acq.spectra - mins ) / ( maxs - mins );      
        end

        % Draw the image to the display
        image( acq.times, f, scaled_spectra ); axis xy;
        drawnow;
 
    else
        % Keep track of skipped displays
        acq.skips = acq.skips + 1;
    end

end

%=======================================================================
% Clean up
%=======================================================================

acq.t1=clock;
elapsed = etime(acq.t1,acq.t0);
fprintf( 'Elapsed time %f seconds\n', elapsed );

% Warn if many skips were encountered
if acq.skips > 5;
    fprintf( '\nWARNING:\nThis program skipped plotting %d times to keep pace.\n', acq.skips )
    fprintf( 'Run again without keyboard interaction or changing the figure size.\n' )
    fprintf( 'If this message reappears you should reduce the plot frequency parameter.\n\n' );
end

% Clean up the analoginput object
if acq.live
  stop(acq.ai);delete( acq.ai );clear acq.ai;
end

% Clean up the figure
delete(fig);
delete(gcf);

if nargout 
   % save up to ten minutes data, preallocated...
    fprintf('Saving output data\n');
   outdata=outdata(1:floor(acq.samples_acquired));
end  

return;

%
%
%=======================================================================
% Functions called
%=======================================================================
%
%


%=======================================================================
% Handle Keypresses
%=======================================================================

% Handle figure window keypress events
function keypress(varargin)

global acq;
keypressed=get(gcf,'CurrentCharacter');

% ignore raw control, shift, alt keys
if keypressed;

    % Save current frame as gif
    if strcmp( keypressed, 'g');
       saveas( acq.fig, sprintf( 'frame%d.png',acq.times(length(acq.times)) ) )
    end
    
    % Offset changes
    increment=1;
    if strcmp( keypressed, 'l');
	acq.offset = acq.offset - increment;
    elseif strcmp( keypressed, 'o');
        acq.offset = acq.offset + increment;

    % Scale changes
    elseif strcmp( keypressed, 'x');
        acq.scale = acq.scale - increment;
    elseif strcmp( keypressed, 's');
        acq.scale = acq.scale + increment;

    % Reset defaults
    elseif strcmp( keypressed, 'd');
        defaults
	acq.restart=1;
    % Normalize spectra
    elseif strcmp( keypressed, 'n');
	   request_normalize
        
    % Quit
    elseif strcmp( keypressed, 'q');
	request_quit

    % Pause
    elseif strcmp( keypressed, 'p');
	   request_pause
  
    % Help
    elseif strcmp( keypressed, 'h');
        audio_instr
  
    % Change colormaps for 0-9,a-c
    elseif strcmp( keypressed, '0' );
        colormap( 'jet' );
    elseif strcmp( keypressed, '1' );
        colormap( 'bone' );
    elseif strcmp( keypressed, '2' );
        colormap( 'colorcube' );
    elseif strcmp( keypressed, '3' );
        colormap( 'cool' );
    elseif strcmp( keypressed, '4' );
        colormap( 'copper' );
    elseif strcmp( keypressed, '5' );
        colormap( 'gray' );
    elseif strcmp( keypressed, '6' );
        colormap( 'hot' );
    elseif strcmp( keypressed, '7' );
        colormap( 'hsv' );
    elseif strcmp( keypressed, '8' );
        colormap( 'autumn' );
    elseif strcmp( keypressed, '9' );
        colormap( 'pink' );
    elseif strcmp( keypressed, 'a' );
        colormap( 'spring' );
    elseif strcmp( keypressed, 'b' );
        colormap( 'summer' );
    elseif strcmp( keypressed, 'c' );
        colormap( 'winter' );        
    end

    update_display
    
end
return

%=======================================================================
% Defaults
%=======================================================================

% Reset defaults
function defaults()
  global acq;
  acq.params.raw_tapers = [2 3];
  acq.moving_window = [0.02 0.02];
  acq.params.tapers=dpsschk(acq.params.raw_tapers,round(acq.params.Fs*acq.moving_window(1)),acq.params.Fs);
  acq.offset = 0;
  acq.scale = 64;
  acq.display_size = 3;
  acq.params.fpass = [50 8000];
  acq.deriv=1;
  acq.log=1;
  acq.bgsub = 1;
  acq.params.pad= 0;
  acq.normalize = 2;
return

function update_display()
global acq;
    set(acq.tapers_ui,'String',sprintf( '%.0f %.0f', acq.params.raw_tapers(1), acq.params.raw_tapers(2) ));
    set(acq.window_ui,'String',sprintf( '%.2f %.2f', acq.moving_window(1), acq.moving_window(2) ));
    set(acq.offset_ui,'String',sprintf( '%d', acq.offset ));
    set(acq.scale_ui,'String',sprintf( '%d', acq.scale ));
    set(acq.display_size_ui,'String',sprintf( '%.1f', acq.display_size ));
    set(acq.frequency_ui,'String',sprintf( '%.1f %.1f', acq.params.fpass(1), acq.params.fpass(2)  ))
    set(acq.derivative_ui,'Value',acq.deriv);
    set(acq.log_ui,'Value',acq.log);
    set(acq.bgsub_ui,'Value',acq.bgsub);
    return


%=======================================================================
% Update ui controls
%=======================================================================

function request_quit(varargin)
	 global acq;
	 acq.stop=1;
return

function request_pause(varargin)
	 global acq;
	 acq.pause = not( acq.pause );
return

function request_normalize(varargin)
	global acq;
        acq.normalize = 2;
return

function update_defaults(varargin)
  global acq;
  defaults
  update_display
  acq.restart=1;
return

function update_tapers(varargin)
	 global acq;
     acq.params.raw_tapers = sscanf(get( gco, 'string' ),'%f %d')';
     acq.params.tapers=dpsschk(acq.params.raw_tapers,round(acq.params.Fs*acq.moving_window(1)),acq.params.Fs); % check tapers
return

function update_window(varargin)
	 global acq;
	 acq.moving_window = sscanf(get( gco, 'string' ),'%f %f');
         acq.params.tapers=dpsschk(acq.params.raw_tapers,round(acq.params.Fs*acq.moving_window(1)),acq.params.Fs);
	 acq.restart = 1;
return

function update_offset(varargin)
	 global acq;
	 acq.offset = sscanf(get( gco, 'string' ),'%f');
	 return

function update_scale(varargin)
	 global acq;
	 acq.scale = sscanf(get( gco, 'string' ),'%f');
	 return

function update_display_size(varargin)
	 global acq;
	 acq.display_size = sscanf(get( gco, 'string' ),'%f');
	 return

function update_fpass(varargin)
	 global acq;
	 acq.params.fpass = sscanf(get( gco, 'string' ),'%f %f');
         acq.restart = 1;
	 return

function update_deriv(varargin)
	 global acq;
     acq.deriv=get( gco, 'Value' );
     acq.normalize=1;
     return
     
function update_log(varargin)
	 global acq;
     acq.log=get( gco, 'Value' );
     acq.normalize=1;
     return

function update_bgsub(varargin)
	 global acq;
     acq.bgsub=get( gco, 'Value' );
     return
     
%=======================================================================
% UI display
%=======================================================================

function fig=create_ui()
    global acq;

    bgcolor = [1 1 1]; % .7 .7 .7
    % ===Create main figure==========================
    fig = figure('Position',centerfig(800,600),...
        'NumberTitle','off',...
        'Name','Real-time spectrogram',...
        'doublebuffer','on',...
        'HandleVisibility','on',...
	'Renderer', 'openGL', ...
	'KeyPressFcn', @keypress, ...
        'Color',bgcolor);

    acq.fig = fig;
    offset = 80;
    % ===text==========
    uicontrol(gcf,'Style','text',...
        'String', 'tapers',...
        'Position',[offset+225 20 45 20],...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'moving win',...
        'Position',[offset+300 20 70 20],...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'offset',...
        'Position',[offset+375 20 30 20],...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'scale',...
        'Position',[offset+410 20 30 20],...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 't axis',...
        'Position',[offset+445 20 30 20],...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'f axis',...
        'Position',[offset+480 20 40 20],...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'deriv',...
        'Position',[offset+550 20 35 20],...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'log',...
        'Position',[offset+580 20 35 20],...
        'BackgroundColor',bgcolor);
    uicontrol(gcf,'Style','text',...
        'String', 'bgsub',...
        'Position',[offset+610 20 35 20],...
        'BackgroundColor',bgcolor);
 
    % ===The quit button===============================
    uicontrol('Style','pushbutton',...
        'Position',[offset+5 5 45 20],...
        'String','Quit',...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'Callback',@request_quit);

    % ===The pause button===============================
    uicontrol('Style','pushbutton',...
        'Position',[offset+55 5 45 20],...
        'String','Pause',...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'Callback',@request_pause);

    % ===The defaults button===============================
    uicontrol('Style','pushbutton',...
        'Position',[offset+105 5 50 20],...
        'String','Defaults',...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'Callback',@update_defaults);

    % ===The normalize button===============================
    uicontrol('Style','pushbutton',...
        'Position',[offset+160 5 60 20],...
        'String','Normalize',...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'Callback',@request_normalize );

    % ===Tapers============================================
    acq.tapers_ui = uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.0f %.0f', acq.params.raw_tapers(1), acq.params.raw_tapers(2) ),...
        'Position',[offset+225 5 70 20],...
        'CallBack', @update_tapers);

    % ===Window============================================
    acq.window_ui=uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.2f %.2f', acq.moving_window(1), acq.moving_window(2) ),...
        'Position',[offset+300 5 70 20],...
        'CallBack', @update_window);

    % ===Offset============================================
    acq.offset_ui = uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%d', acq.offset ),...
        'Position',[offset+375 5 30 20],...
        'CallBack', @update_offset);

    % ===Scale============================================
    acq.scale_ui = uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%d', acq.scale ),...
        'Position',[offset+410 5 30 20],...
        'CallBack', @update_scale);

    % ===display size======================================
    acq.display_size_ui = uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f', acq.display_size ),...
        'Position',[offset+445 5 30 20],...
        'CallBack', @update_display_size);

    % ===frequency axis=====================================
    acq.frequency_ui = uicontrol(gcf,'Style','edit',...
        'String', sprintf( '%.1f %.1f', acq.params.fpass(1), acq.params.fpass(2)  ),...
        'Position',[offset+480 5 80 20],...
        'CallBack', @update_fpass);

    % ===derivative=====================================
    acq.derivative_ui = uicontrol(gcf,'Style','checkbox',...
        'Value',acq.deriv,...
        'Position',[offset+565 5 20 20],...
        'CallBack', @update_deriv);
    
    % ===log=====================================
    acq.log_ui = uicontrol(gcf,'Style','checkbox',...
        'Value',acq.log,...
        'Position',[offset+590 5 20 20],...
        'CallBack', @update_log);

    % ===bgsub=====================================
    acq.bgsub_ui = uicontrol(gcf,'Style','checkbox',...
        'Value',acq.bgsub,...
        'Position',[offset+615 5 20 20],...
        'CallBack', @update_bgsub);

return


%=======================================================================
% Assorted functions
%=======================================================================

function pos = centerfig(width,height)
% Find the screen size in pixels
screen_s = get(0,'ScreenSize');
pos = [screen_s(3)/2 - width/2, screen_s(4)/2 - height/2, width, height];
return


function audio_instr()
% Show instructions

  fprintf('INSTRUCTIONS:\n');
  fprintf('Click on figure window first to activate controls.\n')
  fprintf('Adjust tapers, windows, scales, offsets and axes using the gui\n');
  fprintf('The deriv checkbox toggles derivative of the data\n');
  fprintf('The log checkbox toggles a log of the spectrum\n');
  fprintf('Press d or use defaults button to reset most parameters to defaults.\n')
  fprintf('Press n or use normalize button to normalize spectra based upon values in current display.\n')
  fprintf('Press 0-9,a-c to choose a colormap (default 0).\n')
  fprintf('Press p to pause and unpause display.\n')
  fprintf('Press o and l to adjust offset, or use offset textbox on gui.\n');
  fprintf('Press s and x to adjust scale, or use scale textbox on gui.\n');
  fprintf('Press h for this message.\n')
  fprintf('Press q to quit, or use quit button on gui.\n\n')

return

