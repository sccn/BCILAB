function [patchHandles textHandles] = hlp_vrect(x, varargin)
% Tim Mullen, 2011, SCCN/INC/UCSD

g = finputcheck(varargin, ...
    {'yscale'           'real'      []           1; ...         % vertical scaling (relative to ylimits)
    'dock'              'string'    {'bottom','top','none'} 'top'; ...    % whether to 'dock' the patch at the top or bottom of the axis. If 'none' then scale is location on yaxis of [bottom top] of rect 
    'label'             ''          []              {}; ...
    'axesHandle'        'real'      []              gca; ...
    'textPosition'      'real'      [0 1]           [0.5 0.02]; ...
    'textProperties'    'cell'      []              {}; ...
    'patchProperties'   'cell'      []              {'FaceColor',[0.7 0.7 1],'FaceAlpha',0.5,'EdgeColor',[0.2 0.2 0.2],'EdgeAlpha',0.5} ...
    },'hlp_vrect','error','quiet');

if ischar(g)
    help hlp_vrect;
    error(g);
end

pp = hlp_varargin2struct(g.patchProperties);
if ~isfield(pp,'FaceColor');
    FaceColor = [0.7 0.7 1];
else
    FaceColor = pp.FaceColor;
end
clear pp;

xsize = size(x);
if( xsize(2)~=2 )
    error('x must have two columns');
end

if ~iscell(g.label)
    g.label = {g.label};
end

if length(g.label)~=xsize(1)
    g.label = repmat(g.label,xsize(1),1);
end

if size(g.textPosition,2) ~= 2
    error('textPosition must have two columns [ (x,y) position ]');
end

if size(g.textPosition,1)~=xsize(1)
    g.textPosition=repmat(g.textPosition,xsize(1),1);
end

textPositionX = g.textPosition(:,1);
textPositionY = g.textPosition(:,2);

holdState = ishold(g.axesHandle);
hold(g.axesHandle, 'on');

if strcmpi(g.dock,'none')
    yLimits = g.yscale;
else
    yl = get(g.axesHandle,'ylim');             % Row vector
    yLimits = yl.*g.yscale;            
end

if any(g.yscale<1)
    switch g.dock
        case 'bottom'
            yLimits = yLimits-(yLimits(1)-yl(1))+0.001*diff(yLimits);
        case 'top'
            yLimits = yLimits+(yl(2)-yLimits(2))-0.001*diff(yLimits);
    end
end

% bring desired axes into focus
set(gcf,'CurrentAxes',g.axesHandle);

patchHandles = zeros(1,xsize(1));
textHandles = [];

for k=1:xsize(1)
    % for each patch
    patchHandles(k) = patch([x(k,1) x(k,1) x(k,2) x(k,2)], [yLimits yLimits(end:-1:1)],FaceColor,g.patchProperties{:});
    
    
    if( ~isempty(g.label) )
%         xLimits = get(g.axesHandle,'xlim');
%         xLowerLimit = xLimits(1);
%         xUpperLimit = xLimits(2);
%         xRange      = xUpperLimit - xLowerLimit;
        xPosition   = x(k,1) + textPositionX(k,:)*(x(k,2)-x(k,1));
        
        yUpperLimit = yLimits(2);
        yLowerLimit = yLimits(1);
        yRange      = yUpperLimit - yLowerLimit;
        Yposition   = yLowerLimit + textPositionY(k,:)*yRange;
%         Yposition   = repmat(yPosition, length(x), 1);
        
        % textHandle is a vector correspond to the line
        textHandles(k) = text(xPosition, Yposition, g.label{k}, 'Parent', g.axesHandle,g.textProperties{:});
        % center text
        extent = get(textHandles(k),'Extent');
        set(textHandles(k),'Position',[xPosition-extent(3)/2 Yposition]);
        
    end
end


if( holdState==false )
    hold(g.axesHandle, 'off');
end
