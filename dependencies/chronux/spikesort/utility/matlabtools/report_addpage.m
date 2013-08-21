function report_addpage(filename, fighandle, hangfraction, overlaytext)
%REPORT_ADDPAGE    Builds up a PS report by appending Matlab figures.
%   REPORT_ADDPAGE(FILENAME, FIGUREHANDLE) prints FIGUREHANDLE and appends
%   it to the Postscript file FILENAME (or creates FILENAME if it does not
%   already exist) -- note that the file extension .PS is automatically
%   appended to FILENAME if it is not already there.  If no FIGUREHANDLE
%   is provided (or if FIGUREHANDLE is empty), a blank page is appended to
%   FILENAME.  In either case, the page is printed with a white background
%   on US Letter size paper (8.5-by-11 inches) with 1/4 inch margins.
% 
%   REPORT_ADDPAGE(FILENAME, FIGUREHANDLE, HANGFRACTION) scales the figure
%   such that it occupies HANGFRACTION (range (0.0,1.0]) of a page,
%   spaced such that the whitespace is at the bottom of the page.  This
%   gives the appearance of the end of a section in the report.  If
%   HANGFRACTION is not provided (or is empty), it defaults to 1.0.
%
%   REPORT_ADDPAGE(FILENAME, FIGUREHANDLE, HANGFRACTION, OVERLAYTEXT)
%   also specifies a string caption to be centered and overlaid over the
%   printed page; cell arrays of strings result in multi-line text.  The
%   caption is printed in black, 48-point font centered on the page.  When
%   FIGUREHANDLE is empty, the caption is printed on a blank page.  
%
%   The report uses PostScript rather than the friendlier PDF because of
%   the inflexibility of the Matlab PDF driver (it won't allow pages to be
%   appended).  Use PS2PDF to convert to PDF when the report is complete.
%
%   See also PS2PDF

try    % need a try/catch because we need to cleanup up temp graphics if there's an error

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin  < 2), fighandle = [];  end;
tempfig = 0;  tempaxs = 0;   axshandle = [];
props = {};  oldvals = {};

if (isempty(fighandle)), 
	fighandle = figure;   tempfig = 1;
elseif (~ishandle(fighandle) || ~strcmp(get(fighandle,'type'),'figure'))
	error('Invalid figure handle.');  
end
figure(fighandle);    % make the target the current figure

if (~exist('hangfraction', 'var') || isempty(hangfraction))
	hangfraction = 1.0;
elseif (hangfraction <= 0.0 || hangfraction > 1.0)
	error('The HANGFRACTION must be > 0.0 and <= 1.0.');
end

if (exist('overlaytext', 'var')),    % add overlay text if needed
	tempaxs = 1;
	axshandle = axes('Position', [0 0 1.0 1.0], 'Visible', 'Off');  % 'canvas' covering entire figure
	txthandle = text(0.5, 0.5, overlaytext, 'FontSize', 48, 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Prep the Figure %%%%%%%%%%%%%%%%%%%%%%%%%%
props = {'PaperUnits', 'PaperSize', 'PaperType', 'InvertHardCopy', 'Color'};
newvals = {'inches', [8.5 11], 'usletter', 'off', 'w'};
oldvals = get(fighandle, props);

propvals = cat(1, props, newvals);
set(fighandle, propvals{:});

set(gcf, 'PaperPosition', [0.25 10.5*(1-hangfraction) 8 10.5*hangfraction]);


%%%%%%%%%%%%%%%%%%%%%%%%%% Print the Figure %%%%%%%%%%%%%%%%%%%%%%%%%%
print('-dpsc2', '-append', '-painters', filename);  % painters renderer is good for EPS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clean Up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cleanup(tempfig, fighandle, tempaxs, axshandle, props, oldvals);

catch
	try   cleanup(tempfig, fighandle, tempaxs, axshandle, props, oldvals);
    catch 
	end	
	rethrow(lasterror);
end


% Cleans up any temporary objects so people get their graphics in good shape.
function cleanup(tempfig, fighandle, tempaxs, axshandle, props, oldvals)

if (tempaxs && ishandle(axshandle)), delete(axshandle);  end;  % deletes caption too, if it exists
if (tempfig && ishandle(fighandle))
	close(fighandle);
else
	propvals = cat(1, props, oldvals);
	set(fighandle, propvals{:});
end
