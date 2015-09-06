function F = stoploop(str)
% STOPLOOP - creates stop button to have a user interrupt a loop
% 
%  FS = STOPLOOP creates a message box window and returns a structure FS that
%  holds two functions, called FS.Stop and FS.Clear. The function FS.Stop()
%  will return true, if the OK button has been clicked (or the message box
%  has been removed), so that a loop can be interrupted.
%  The function FS.Clear() can be used to remove the message box, if a loop
%  has ended without user interruption.
%
%  FS = STOPLOOP(STR) uses the string STR to display instead of the default
%  'Stop the Loop'.
%
%  Example:
%    tic ;         % We will measure elapsed time in a loop
%                  % Set up the stop box:
%    FS = stoploop({'Stop me before', '5 seconds have elapsed'}) ;
%                  % Display elapsed time
%    fprintf('\nSTOPLOOP: elapsed time (s): %5.2f\n',toc)
%                  % start the loop
%    while(~FS.Stop() && toc < 5),       % Check if the loop has to be stopped
%        fprintf('%c',repmat(8,6,1)) ;   % clear up previous time
%        fprintf('%5.2f\n',toc) ;        % display elapsed time
%    end
%    FS.Clear() ;  % Clear up the box
%    clear FS ;    % this structure has no use anymore
%
%  Notes:
%  - The function call F.Stop() issues a drawnow command.
%  - The function call F.Clear() does nothing when the message box already
%    has been cleared.
%  
%  See also MSGBOX, KEYBOARD, DBSTOP, WHILE, FOR
%           GETKEYNOW, GETKEYWAIT, GETKEY (on the File Exchange)

% for Matlab R14+
% version 1.0 (jun 2008)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History:
% 1.0 (jun 2008) Inspired by several requests over the years on CSSM

error(nargoutchk(1,1,nargout)) ;

if nargin,
    if ~ischar(str) && ~iscellstr(str),
        error([mfilename ':InputString'],...
            'Input should be a string, or a cell array of strings.') ;
    end
else
    % default message string
    str = 'Stop the Loop' ;
end

% create a msgbox displaying the string
H = msgbox(str,'STOPLOOP') ;

% create the two anonymous functions
F.Stop = @() stopfun(H) ; % false if message box still exists
F.Clear = @() clearfun(H) ; % delete message box

function r = stopfun(H)
drawnow ;          % ensure that button presses are recorded
r = ~ishandle(H) ; % false if message box still exists

function clearfun(H)
% clear the message box if it still exists
if ishandle(H),
    delete(H) ;
end


