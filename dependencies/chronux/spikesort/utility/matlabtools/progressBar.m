function progressBar(progressFraction, updateEvery, infoString)
%PROGRESSBAR       Improved progress indicator.
%   PROGRESSBAR(FRACTION, UPDATEEVERY, INFOSTRING) modifies Matlab's
%   waitbar functionality to allow the indicator to decrease, automate
%   creation/deletion and restrict updates so the indicator is only
%   redrawn every kth iteration.
%
%   If the indicator figure does not exist when PROGRESSBAR is called, it
%   will be created (via WAITBAR) initialized with the message INFOSTRING
%   and percent complete given by FRACTION * 100%.
%
%   If the figure already exists, it is updated if FRACTION < 1.0, and if
%   the number of calls since initialization is an integral multiple of
%   the UPDATEEVERY parameter used to initialize the figure (this update
%   rate can be changed by passing in a new value).  Otherwise, nothing
%   happens.
%
%   To force an update on an existing progressbar regardless of the the
%   UPDATEEVERY count, pass in [] for the second argument.  (E.g., this
%   can be used to force info string updates.)
%   
%   Finally, if FRACTION is >= 1.0, the indicator is deleted (or it is not
%   shown if it does not already exist).
%
%   See also WAITBAR.

persistent progressHandle iter updateRate;

if (progressFraction >= 1.0)  % delete if it exists & reset the handle
    if (~isempty(progressHandle) && ishandle(progressHandle))
        delete(progressHandle);
    end
    clear progressHandle iter;
elseif (isempty(progressHandle) || ~ishandle(progressHandle))  % create if it doesn't exist
    if (nargin < 3),  infoString = '';  end;
    if (nargin < 2)
        error('At least two arguments are required when initializing.');
    end
    updateRate = updateEvery;
    iter = 1;
    progressHandle = waitbar(progressFraction, infoString);
	posDef = get(0, 'DefaultFigurePosition');  posWait = get(progressHandle, 'Position');
	set(progressHandle, 'Position', [posDef(1:2) posWait(3:4)]);
	
    p = findobj(progressHandle, 'Type', 'Patch');
    set(p, 'EraseMode', 'normal');   % allows the indicator to shrink properly    
else  % update it if thats the right thing to do
    if (nargin > 1)
        if (isempty(updateRate))
            update = 1;
        else
            updateRate = updateEvery;
            update = (iter == updateRate);
        end
    else
        update = (iter == updateRate);
    end
    if (update)
        iter = 1;
        if (nargin < 3)
            waitbar(progressFraction, progressHandle);
        else
            waitbar(progressFraction, progressHandle, infoString);
        end
    else
        iter = iter + 1;
    end    
end

