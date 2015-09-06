function cancel = multiWaitbar( label, varargin )
%multiWaitbar: add, remove or update an entry on the multi waitbar
%
%   multiWaitbar(LABEL,VALUE) adds a waitbar for the specified label, or
%   if it already exists updates the value. LABEL must be a string and
%   VALUE a number between zero and one or the string 'Close' to remove the
%   entry Setting value equal to 0 or 'Reset' will cause the progress bar
%   to reset and the time estimate to be re-initialised.
%
%   multiWaitbar(LABEL,COMMAND,VALUE,...) passes one or more command/value
%   pairs for changing the named waitbar entry. Possible commands include:
%   'Value'       Set the value of the named waitbar entry. The
%                 corresponding value must be a number between 0 and 1.
%   'Increment'   Increment the value of the named waitbar entry. The
%                 corresponding value must be a number between 0 and 1.
%   'Color'       Change the color of the named waitbar entry. The
%                 value must be an RGB triple, e.g. [0.1 0.2 0.3].
%   'Reset'       Set the named waitbar entry back to zero and reset its
%                 timer. No value need be specified.
%   'CanCancel'   [on|off] should a "cancel" button be shown for this bar
%                 (default 'off').
%   'CancelFcn'   Function to call in the event that the user cancels.
%   'ResetCancel' Reset the "cancelled" flag for an entry (ie. if you
%                 decide not to cancel).
%   'Close'       Remove the named waitbar entry.
%
%   cancel = multiWaitbar(LABEL,VALUE) also returns whether the user has
%   clicked the "cancel" button for this entry (true or false). Two
%   mechanisms are provided for cancelling an entry if the 'CanCancel'
%   setting is 'on'. The first is just to check the return argument and if
%   it is true abort the task. The second is to set a 'CancelFcn' that is
%   called when the user clicks the cancel button, much as is done for
%   MATLAB's built-in WAITBAR. In either case, you can use the
%   'ResetCancel' command if you don't want to cancel afterall. 
%
%   multiWaitbar('CLOSEALL') closes the waitbar window.
%
%   Example:
%   multiWaitbar( 'CloseAll' );
%   multiWaitbar( 'Task 1', 0 );
%   multiWaitbar( 'Task 2', 0.5, 'Color', [0.2 0.6 0.2] );
%   multiWaitbar( 'Task 1', 'Value', 0.1 );
%   multiWaitbar( 'Task 2', 'Increment', 0.2 );
%   multiWaitbar( 'Task 2', 'Close' );
%   multiWaitbar( 'Task 1', 'Reset' );
%
%   Example:
%   multiWaitbar( 'Task 1', 0, 'CancelFcn', @(a,b) disp( ['Cancel ',a] ) );
%   for ii=1:100
%      abort = multiWaitbar( 'Task 1', ii/100 );
%      if abort
%         % Here we would normally ask the user if they're sure
%         break
%      else
%         pause( 1 )
%      end
%   end
%   multiWaitbar( 'Task 1', 'Close' )
%
%   Example:
%   multiWaitbar( 'CloseAll' );
%   multiWaitbar( 'Red...',    7/7, 'Color', [0.8 0.0 0.1] );
%   multiWaitbar( 'Orange...', 6/7, 'Color', [1.0 0.4 0.0] );
%   multiWaitbar( 'Yellow...', 5/7, 'Color', [0.9 0.8 0.2] );
%   multiWaitbar( 'Green...',  4/7, 'Color', [0.2 0.9 0.3] );
%   multiWaitbar( 'Blue...',   3/7, 'Color', [0.1 0.5 0.8] );
%   multiWaitbar( 'Indigo...', 2/7, 'Color', [0.4 0.1 0.5] );
%   multiWaitbar( 'Violet...', 1/7, 'Color', [0.8 0.4 0.9] );

%   Author: Ben Tordoff
%   Copyright 2007-2010 The MathWorks Ltd
%   $Revision: 47$
%   $Date: 2007-10-17$

persistent figh;

% Check basic inputs
error( nargchk( 1, inf, nargin ) );
if ~ischar( label )
    error( 'multiWaitbar:BadArg', 'LABEL must be the name of the progress entry (i.e. a string)' );
end

% Try to get hold of the figure
if isempty( figh ) || ~ishandle( figh )
    figh = findall( 0, 'Type', 'figure', 'Tag', 'multiWaitbar:Figure' );
    if isempty(figh)
        figh = iCreateFig();
    else
        figh = handle( figh(1) );
    end
end

% Check for close all and stop early
if any( strcmpi( label, {'CLOSEALL','CLOSE ALL'} ) )
    delete( figh );
    return;
end

% Make sure we're ons-screen
if ~strcmpi( figh.Visible, 'on' )
    figh.Visible = 'on';
end

% Get the list of entries and see if this one already exists
entries = getappdata( figh, 'ProgressEntries' );
if isempty(entries)
    idx = [];
else
    idx = find( strcmp( label, {entries.Label} ), 1, 'first' );
end
bgcol = getappdata( figh, 'DefaultProgressBarBackgroundColor' );

% If it doesn't exist, create it
needs_redraw = false;
entry_added = isempty(idx);
if entry_added
    % Create a new entry
    defbarcolor = getappdata( figh, 'DefaultProgressBarColor' );
    entries = iAddEntry( figh, entries, label, 0, defbarcolor, bgcol );
    idx = numel( entries );
end

% Check if the user requested a cancel
if nargout
    cancel = entries(idx).Cancel;
end

% Parse the inputs. We shortcut the most common case as an efficiency
force_update = false;
if nargin==2 && isnumeric( varargin{1} )
    entries(idx).LastValue = entries(idx).Value;
    entries(idx).Value = max( 0, min( 1, varargin{1} ) );
    needs_update = true;
else
    [params,values] = iParseInputs( varargin{:} );
    
    needs_update = false;
    for ii=1:numel( params )
        switch upper( params{ii} )
            case 'VALUE'
                entries(idx).LastValue = entries(idx).Value;
                entries(idx).Value = max( 0, min( 1, values{ii} ) );
                needs_update = true;
                
            case {'INC','INCREMENT'}
                entries(idx).LastValue = entries(idx).Value;
                entries(idx).Value = max( 0, min( 1, entries(idx).Value + values{ii} ) );
                needs_update = true;
                
            case {'COLOR','COLOUR'}
                entries(idx).CData = iMakeColors( values{ii}, 16 );
                needs_update = true;
                force_update = true;
                
            case {'CANCANCEL'}
                if ~ischar( values{ii} ) || ~any( strcmpi( values{ii}, {'on','off'} ) )
                    error( 'multiWaitbar:BadString', 'Parameter ''CanCancel'' must be a ''on'' or ''off''.' );
                end
                entries(idx).CanCancel = strcmpi( values{ii}, 'on' );
                needs_redraw = true;
                
            case {'CANCELFCN'}
                if ~isa( values{ii}, 'function_handle' )
                    error( 'multiWaitbar:BadFunction', 'Parameter ''CancelFcn'' must be a valid function handle.' );
                end
                entries(idx).CancelFcn = values{ii};
                if ~entries(idx).CanCancel
                    entries(idx).CanCancel = true;
                end
                needs_redraw = true;
                
            case {'CLOSE','DONE'}
                if ~isempty(idx)
                    % Remove the selected entry
                    entries = iDeleteEntry( entries, idx );
                end
                if isempty( entries )
                    delete( figh );
                    % With the window closed, there's nothing else to do
                    return;
                else
                    needs_redraw = true;
                end
                % We can't continue after clearing the entry, so jump out
                break;
            case {'RESETCANCEL'}
                % Set the cancel flag so that the user is told on next update
                entries(idx).Cancel = false;
                setappdata( figh, 'ProgressEntries', entries );
                
            otherwise
                error( 'multiWaitbar:BadArg', 'Unrecognised command: ''%s''', params{ii} );
                
        end
    end
end
if needs_redraw
    setappdata( figh, 'ProgressEntries', entries );
    iRedraw( figh );
    % NB: Redraw includes updating all bars, so never need to do both
elseif needs_update
    [entries(idx),updated] = iUpdateEntry( entries(idx), force_update );
    setappdata( figh, 'ProgressEntries', entries );
    if updated && ~entry_added
        drawnow('expose');
    end
end
if entry_added || needs_redraw
    % If the shape or size has changed, do a full redraw, including events
    drawnow();
end

end % multiWaitbar


%-------------------------------------------------------------------------%
function [params, values] = iParseInputs( varargin )
% Parse the input arguments, extracting a list of commands and values
idx = 1;
params = {};
values = {};
if nargin==0
    return;
end
if isnumeric( varargin{1} )
    params{idx} = 'Value';
    values{idx} = varargin{1};
    idx = idx + 1;
end

while idx <= nargin
    param = varargin{idx};
    if ~ischar( param )
        error( 'multiWaitbar:BadSyntax', 'Additional properties must be supplied as property-value pairs' );
    end
    params{end+1,1} = param; %#ok<AGROW>
    values{end+1,1} = []; %#ok<AGROW>
    switch upper( param )
        case {'DONE','CLOSE'}
            % No value needed, and stop
            break;
        case {'RESET','ZERO','SHOW'}
            params{end} = 'Value'; %#ok<AGROW>
            values{end} = 0; %#ok<AGROW>
            idx = idx + 1;
        otherwise
            if idx==nargin
                error( 'multiWaitbar:BadSyntax', 'Additional properties must be supplied as property-value pairs' );
            end
            values{end,1} = varargin{idx+1}; %#ok<AGROW>
            idx = idx + 2;
    end
end
if isempty( params )
    error( 'multiWaitbar:BadSyntax', 'Must specify a value or a command' );
end
end % iParseInputs

%-------------------------------------------------------------------------%
function fobj = iCreateFig()
try
    icadefs;
catch err
	GUIBACKCOLOR        =   get(0,'DefaultUIControlBackgroundColor');
    GUITEXTCOLOR        =   [0 0 0];
end
% Create the progress bar group window
bgcol = GUIBACKCOLOR; %get(0,'DefaultUIControlBackgroundColor');
f = figure( ...
    'Name', 'Progress', ...
    'Tag', 'multiWaitbar:Figure', ...
    'Color', bgcol, ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'HandleVisibility', 'off', ...
    'IntegerHandle', 'off', ...
    'Visible', 'off', ...
    'NumberTitle', 'off' );
% Resize
fobj = handle( f );
fobj.Position(3:4) = [360 42];
setappdata( fobj, 'ProgressEntries', [] );
% Make sure we have the image
defbarcolor = [0.8 0.0 0.1];
barbgcol = uint8( 255*0.75*bgcol );
setappdata( fobj, 'DefaultBackgroundColor', GUIBACKCOLOR );
setappdata( fobj, 'DefaultTextColor', GUITEXTCOLOR );
setappdata( fobj, 'DefaultProgressBarBackgroundColor', barbgcol );
setappdata( fobj, 'DefaultProgressBarColor', defbarcolor );
setappdata( fobj, 'DefaultProgressBarSize', [350 16] );
% Setup the resize function after we've finished setting up the figure to
% avoid excessive redraws
fobj.ResizeFcn = @iRedraw;
fobj.CloseRequestFcn = @iCloseFigure;
end % iCreateFig

%-------------------------------------------------------------------------%
function cdata = iMakeColors( baseColor, height )
% Creates a shiny bar from a single base color
lightColor = [1 1 1];

if isa( baseColor, 'uint8' )
    baseColor = double( baseColor ) / 255;
end
cols = repmat( baseColor, [height, 1] );

breaks = max( 1, round( height * [1 25 50 75 88 100] / 100 ) );
cols(breaks(1),:) = 0.6*baseColor;
cols(breaks(2),:) = lightColor - 0.4*(lightColor-baseColor);
cols(breaks(3),:) = baseColor;
cols(breaks(4),:) = min( baseColor*1.2, 1.0 );
cols(breaks(5),:) = min( baseColor*1.4, 0.95 ) + 0.05;
cols(breaks(6),:) = min( baseColor*1.6, 0.9 ) + 0.1;

y = 1:height;
cols(:,1) = max( 0, min( 1, interp1( breaks, cols(breaks,1), y, 'cubic' ) ) );
cols(:,2) = max( 0, min( 1, interp1( breaks, cols(breaks,2), y, 'cubic' ) ) );
cols(:,3) = max( 0, min( 1, interp1( breaks, cols(breaks,3), y, 'cubic' ) ) );
cdata = uint8( 255 * cat( 3, cols(:,1), cols(:,2), cols(:,3) ) );
end % iMakeColors


%-------------------------------------------------------------------------%
function cdata = iMakeBackground( baseColor, height )
% Creates a shaded background
if isa( baseColor, 'uint8' )
    baseColor = double( baseColor ) / 255;
end

ratio = 1 - exp( -0.5-2*(1:height)/height )';
cdata = uint8( 255 * cat( 3, baseColor(1)*ratio, baseColor(2)*ratio, baseColor(3)*ratio ) );
end % iMakeBackground

%-------------------------------------------------------------------------%
function entries = iAddEntry( parent, entries, label, value, color, bgcolor )
% Add a new entry to the progress bar

% Create bar coloring
psize = getappdata( parent, 'DefaultProgressBarSize' );
cdata = iMakeColors( color, 16 );
% Create background image
barcdata = iMakeBackground( bgcolor, psize(2) );

% Work out the size in advance
labeltext = uicontrol( 'Style', 'Text', ...
    'String', label, ...
    'Parent', parent, ...
    'HorizontalAlignment', 'Left' , ...
    'BackgroundColor', getappdata( parent, 'DefaultBackgroundColor' ), ...
    'ForegroundColor', getappdata( parent, 'DefaultTextColor' ));
etatext = uicontrol( 'Style', 'Text', ...
    'String', '', ...
    'Parent', parent, ...
    'HorizontalAlignment', 'Right' , ...
    'BackgroundColor', getappdata( parent, 'DefaultBackgroundColor' ), ...
    'ForegroundColor', getappdata( parent, 'DefaultTextColor' ));
progresswidget = uicontrol( 'Style', 'Checkbox', ...
    'String', '', ...
    'Parent', parent, ...
    'Position', [5 5 psize], ...
    'CData', barcdata , ...
    'BackgroundColor', getappdata( parent, 'DefaultBackgroundColor' ), ...
    'ForegroundColor', getappdata( parent, 'DefaultTextColor' ));
cancelwidget = uicontrol( 'Style', 'PushButton', ...
    'String', '', ...
    'FontWeight', 'Bold', ...
    'Parent', parent, ...
    'Position', [5 5 16 16], ...
    'CData', iMakeCross( 8 ), ...
    'Callback', @(src,evt) iCancelEntry( src, label ), ...
    'Visible', 'off' );
mypanel = uipanel( 'Parent', parent, 'Units', 'Pixels' , ...
    'BackgroundColor', getappdata( parent, 'DefaultBackgroundColor' ), ...
    'ForegroundColor', getappdata( parent, 'DefaultTextColor' ));

newentry = struct( ...
    'Label', label, ...
    'Value', value, ...
    'LastValue', inf, ...
    'Created', tic(), ...
    'LabelText', labeltext, ...
    'ETAText', etatext, ...
    'ETAString', '', ...
    'Progress', progresswidget, ...
    'ProgressSize', psize, ...
    'Panel', mypanel, ...
    'BarCData', barcdata, ...
    'CData', cdata, ...
    'BackgroundCData', barcdata, ...
    'CanCancel', false, ...
    'CancelFcn', [], ...
    'CancelButton', cancelwidget, ...
    'Cancel', false );
if isempty( entries )
    entries = newentry;
else
    entries = [entries;newentry];
end
% Store in figure before the redraw
setappdata( parent, 'ProgressEntries', entries );
if strcmpi( get( parent, 'Visible' ), 'on' )
    iRedraw( parent, [] );
else
    set( parent, 'Visible', 'on' );
end
end % iAddEntry

%-------------------------------------------------------------------------%
function entries = iDeleteEntry( entries, idx )
delete( entries(idx).LabelText );
delete( entries(idx).ETAText );
delete( entries(idx).CancelButton );
delete( entries(idx).Progress );
delete( entries(idx).Panel );
entries(idx,:) = [];
end % iDeleteEntry

%-------------------------------------------------------------------------%
function entries = iCancelEntry( src, name )
figh = ancestor( src, 'figure' );
entries = getappdata( figh, 'ProgressEntries' );
if isempty(entries)
    % The entries have been lost - nothing can be done.
    return
end
idx = find( strcmp( name, {entries.Label} ), 1, 'first' );

% Set the cancel flag so that the user is told on next update
entries(idx).Cancel = true;
setappdata( figh, 'ProgressEntries', entries );

% If a user function is supplied, call it
if ~isempty( entries(idx).CancelFcn )
    feval( entries(idx).CancelFcn, name, 'Cancelled' );
end

end % iCancelEntry


%-------------------------------------------------------------------------%
function [entry,updated] = iUpdateEntry( entry, force )
% Update one progress bar

% Some constants
marker_weight = 0.8;
shadow1_weight = 0.4;
shadow2_weight = 0.7;

% Check if the label needs updating
updated = force;
val = entry.Value;
lastval = entry.LastValue;

% Now update the bar
psize = entry.ProgressSize;
filled = max( 1, round( val*psize(1) ) );
lastfilled = max( 1, round( lastval*psize(1) ) );

% We do some careful checking so that we only redraw what we have to. This
% makes a small speed difference, but every little helps!
if force || (filled<lastfilled)
    % Create the bar background
    bgim = entry.BackgroundCData(:,ones( 1, psize(1)-filled ),:);
    % We use slightly bizarre indexing notation to achieve REPMAT of the
    % column at roughly 10x the speed. Blame Jon Cherrie (who showed me
    % that it was even faster than BSXFUN at doing REPMAT)!
    barim = entry.CData(:,ones( 1, filled ),:);
    progresscdata = [barim,bgim];
    
    % Add light/shadow around the markers
    markers = round( (0.1:0.1:val)*psize(1) );
    highlight = [marker_weight*entry.CData, 255 - marker_weight*(255-entry.CData)];
    for ii=1:numel( markers )
        progresscdata(:,markers(ii)+[-1,0],:) = highlight;
    end
    
    % Add highlight and shadow to the ends of the bar
    progresscdata(:,1,:) = 255 - shadow1_weight*(255-entry.CData);
    progresscdata(:,filled,:) = shadow1_weight*entry.CData;
    if filled>=4
        progresscdata(:,2,:) = 255 - shadow2_weight*(255-entry.CData);
        progresscdata(:,filled-1,:) = shadow2_weight*entry.CData;
    end
    
    % Set the image into the checkbox
    entry.BarCData = progresscdata;
    set( entry.Progress, 'cdata', progresscdata );
    updated = true;
    
elseif filled > lastfilled
    % Just need to update the existing data
    progresscdata = entry.BarCData;
    startidx = max(1,lastfilled-1);
    % Repmat is the obvious way to fill the bar, but BSXFUN is often
    % faster. Indexing is obscure but faster still.
    progresscdata(:,startidx:filled,:) = entry.CData(:,ones(1,filled-startidx+1),:);
    
    % Add light/shadow around the markers
    markers = round( (0.1:0.1:val)*psize(1) );
    markers(markers<startidx) = [];
    highlight = [marker_weight*entry.CData, 255 - marker_weight*(255-entry.CData)];
    for ii=1:numel( markers )
        progresscdata(:,markers(ii)+[-1,0],:) = highlight;
    end
    
    % Add highlight and shadow to the ends of the bar
    if filled>2 && lastfilled<=3
        progresscdata(:,2,:) = 255 - shadow2_weight*(255-entry.CData);
        progresscdata(:,1,:) = 255 - shadow1_weight*(255-entry.CData);
    end
    progresscdata(:,filled,:) = shadow1_weight*entry.CData;
    if filled>=4
        progresscdata(:,filled-1,:) = shadow2_weight*entry.CData;
    end
    entry.BarCData = progresscdata;
    set( entry.Progress, 'CData', progresscdata );
    updated = true;
end

% Now work out the remaining time
minTime = 3; % secs
if val <= 0
    % Zero value, so clear the eta
    entry.Created = tic();
    elapsedtime = 0;
    eta = '';
else
    elapsedtime = toc( entry.Created ); % in seconds
    
    % Only show the remaining time if we've had time to estimate
    if elapsedtime < minTime
        % Not enough time has passed since starting, so leave blank
        eta = '';
    else
        % Calculate a rough ETA
        remainingtime = elapsedtime * (1-val) / val;
        
        if remainingtime > 172800 % 2 days
            eta = sprintf( '%d days', round(remainingtime/86400) );
        else
            if remainingtime > 7200 % 2 hours
                eta = sprintf( '%d hours', round(remainingtime/3600) );
            else
                if remainingtime > 120 % 2 mins
                    eta = sprintf( '%d mins', round(remainingtime/60) );
                else
                    % Seconds
                    remainingtime = round( remainingtime );
                    if remainingtime > 1
                        eta = sprintf( '%d secs', remainingtime );
                    elseif remainingtime == 1
                        eta = '1 sec';
                    else
                        eta = ''; % Nearly done (<1sec)
                    end
                end
            end
        end
    end
end

if ~isequal( eta, entry.ETAString )
    set( entry.ETAText, 'String', eta );
    entry.ETAString = eta;
    updated = true;
end

% Update the label too
if elapsedtime > minTime
    decval = round( val*100 );
    if force || (decval ~= round( lastval*100 ))
        labelstr = [entry.Label, sprintf( ' (%d%%)', decval )];
        set( entry.LabelText, 'String', labelstr );
        updated = true;
    end
end

end % iUpdateEntry

%-------------------------------------------------------------------------%
function iCloseFigure( fig, evt ) %#ok<INUSD>
% Closing the figure just makes it invisible
set( fig, 'Visible', 'off' );
end % iCloseFigure

%-------------------------------------------------------------------------%
function iRedraw( fig, evt ) %#ok<INUSD>
entries = getappdata( fig, 'ProgressEntries' );
fobj = handle( fig );
p = fobj.Position;
% p = get( fig, 'Position' );
border = 5;
textheight = 16;
barheight = 16;
panelheight = 10;
N = max( 1, numel( entries ) );

% Check the height is correct
heightperentry = textheight+barheight+panelheight;
requiredheight = 2*border + N*heightperentry - panelheight;
if ~isequal( p(4), requiredheight )
    p(2) = p(2) + p(4) - requiredheight;
    p(4) = requiredheight;
    % Changing the position will re-fire this callback, so return
    % immediately
    set( fig, 'Position', p )
    return;
end
ypos = p(4) - border;
width = p(3) - 2*border;
setappdata( fig, 'DefaultProgressBarSize', [width barheight] );

for ii=1:numel( entries )
    set( entries(ii).LabelText, 'Position', [border ypos-textheight width*0.75 textheight] );
    set( entries(ii).ETAText, 'Position', [border+width*0.75 ypos-textheight width*0.25 textheight] );
    ypos = ypos - textheight;
    if entries(ii).CanCancel
        set( entries(ii).Progress, 'Position', [border ypos-barheight width-barheight+1 barheight] );
        entries(ii).ProgressSize = [width-barheight barheight];
        set( entries(ii).CancelButton, 'Visible', 'on', 'Position', [p(3)-border-barheight ypos-barheight barheight barheight] );
    else
        set( entries(ii).Progress, 'Position', [border ypos-barheight width+1 barheight] );
        entries(ii).ProgressSize = [width barheight];
        set( entries(ii).CancelButton, 'Visible', 'off' );
    end
    ypos = ypos - barheight;
    set( entries(ii).Panel, 'Position', [-500 ypos-500-panelheight/2 p(3)+1000 500] );
    ypos = ypos - panelheight;
    entries(ii) = iUpdateEntry( entries(ii), true );
end
setappdata( fig, 'ProgressEntries', entries );
end % iRedraw

function cdata = iMakeCross( sz )
% Create a cross-shape image

cdata = nan( sz, sz );
for ii=1:sz
    cdata(ii,ii) = 0;
    cdata(sz-ii+1,ii) = 0;
end
for ii=2:sz
    cdata(ii,ii-1) = 0;
    cdata(ii-1,ii) = 0;
    cdata(sz-ii+1,ii-1) = 0;
    cdata(ii,sz-ii+2) = 0;
end
cdata = cat( 3, cdata, cdata, cdata );

end % iMakeCross
    
