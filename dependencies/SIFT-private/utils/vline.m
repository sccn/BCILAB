function [lineHandles textHandles] = vline(x, lineType, label, textPosition, axesHandle, textColor)

% Draws vectical lines (specified by a vector of x)
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
% Acknowledgement:  It was based on vline() written by Brandon Kuczenski.

    if( ~isvector(x) )
        error('x must be a vector');
    else
        x=x(:);     % Make sure it's column vector
    end

    if( ~exist('label', 'var') )        label = [];         end
    if( ~exist('lineType', 'var') )     lineType = [];      end
    if( ~exist('axesHandle', 'var') )   axesHandle = gca;   end
    if( ~exist('textPosition', 'var') ) textPosition = [];  end
    if( ~exist('textColor', 'var') )    textColor = [];     end
    
    if( isempty(axesHandle) )           axesHandle = gca;   end
    
    if( isempty(textPosition) )
        textPositionX = 0.02;
        textPositionY = 0.02;
    elseif( isscalar(textPosition) )
        textPositionX = 0.02;
        textPositionY = textPosition; 
    elseif( length(textPosition)~=2 )
        error('Invalid textPosition');
    else
        textPositionX = textPosition(1);
        textPositionY = textPosition(2);
    end
    clear textPosition;     % Avoid typos for similarly named variables
    
    holdState = ishold(axesHandle);
    hold(axesHandle, 'on');
    
    yLimits=get(axesHandle,'ylim');             % Row vector    
    Ylimits=repmat(yLimits', 1, length(x));
    
    % Example: for horizontal lines
    % X = [2 5        Y = [3 3
    %      2 5];           4 4];
            
    if( isempty(lineType) )
        lineHandles = plot(axesHandle, [x';x'], Ylimits);
    else
        lineHandles = plot(axesHandle, [x';x'], Ylimits, lineType);
    end
    
    if( ~isempty(label) )
        xLimits = get(axesHandle,'xlim');
        xLowerLimit = xLimits(2);
        xUpperLimit = xLimits(1);        
        xRange      = xUpperLimit - xLowerLimit;
        xPosition   = x+textPositionX*xRange;
        
        yUpperLimit = yLimits(2);
        yLowerLimit = yLimits(1);
        yRange      = yUpperLimit - yLowerLimit;
        yPosition   = yLowerLimit+textPositionY*yRange;
        Yposition   = repmat(yPosition, length(x), 1);
        
        % textHandle is a vector correspond to the line
        textHandles = text(xPosition, Yposition, label, 'Parent', axesHandle);
    
        for k=1:length(x)
            % Set the text colors to be identical to line colors
            if isempty(textColor)
                set( textHandles(k), 'color', get(lineHandles(k), 'color') );
            else
                set( textHandles(k), 'color', textColor );
            end
        end                            
    end
    
    if( holdState==false )
        hold(axesHandle, 'off');
    end
    % this last part is so that it doesn't show up on legends
    set(lineHandles,'tag','vline','handlevisibility','off')
