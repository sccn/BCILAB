function [lineHandles textHandles] = hline(y, lineType, label, textPosition, axesHandle)

% Draws horizontal lines (specified by a vector of y)
%
% Except the first argument 'y', everything else are optional. 
% To preserve input argument order, enter [] for unused input arguments
%
% lineType:     same as the ones used in plot(...)
% label:        simple strings applies to all lines
%               cell strings applies to each line
% textPosition: [x, y] bewteen 0 and 1 relative to the 'lower' end
%
% Author:      Hoi Wong (wonghoi.ee@gmail.com)
% Date:        02/14/2008
%
% Acknowledgement:  It was based on hline() written by Brandon Kuczenski.

    if( ~isvector(y) )
        error('y must be a vector');
    else
        y=y(:);     % Make sure it's column vector
    end

    if( ~exist('label', 'var') )        label = [];         end
    if( ~exist('lineType', 'var') )     lineType = [];      end
    if( ~exist('axesHandle', 'var') )   axesHandle = gca;   end
    if( ~exist('textPosition', 'var') ) textPosition = [];  end
    
    if( isempty(axesHandle) )           axesHandle = gca;   end
    
    if( isempty(textPosition) )
        textPositionX = 0.02;
        textPositionY = 0.02;
    elseif( isscalar(textPosition) )
        textPositionX = textPosition;
        textPositionY = 0.02; 
    elseif( length(textPosition)~=2 )
        error('Invalid textPosition');
    else
        textPositionX = textPosition(1);
        textPositionY = textPosition(2);
    end
    clear textPosition;     % Avoid typos for similarly named variables
               
    holdState = ishold(axesHandle);
    hold(axesHandle, 'on');
    
    xLimits=get(axesHandle,'xlim');             % Row vector    
    Xlimits=repmat(xLimits', 1, length(y));
    
    % Example: for horizontal lines
    % X = [2 2        Y = [3 4
    %      5 5];           3 4];
            
    if( isempty(lineType) )
        lineHandles = plot(axesHandle, Xlimits, [y';y']);
    else
        lineHandles = plot(axesHandle, Xlimits, [y';y'], lineType);
    end
    
    if( ~isempty(label) )
        yLimits = get(axesHandle,'ylim');
        yLowerLimit = yLimits(2);
        yUpperLimit = yLimits(1);        
        yRange      = yUpperLimit - yLowerLimit;
        yPosition   = y-textPositionY*yRange;
        
        xUpperLimit = xLimits(2);
        xLowerLimit = xLimits(1);
        xRange      = xUpperLimit - xLowerLimit;
        xPosition   = xLowerLimit+textPositionX*xRange;
        Xposition   = repmat(xPosition, length(y), 1);
        
        % textHandle is a vector correspond to the line
        textHandles = text(Xposition, yPosition, label, 'Parent', axesHandle);
        
        for k=1:length(y)
            % Set the text colors to be identical to line colors
            set( textHandles(k), 'color', get(lineHandles(k), 'color') );
        end                    
    end
    
    if( holdState==false )
        hold(axesHandle, 'off');
    end
    % this last part is so that it doesn't show up on legends
    set(lineHandles,'tag','hline','handlevisibility','off') 
