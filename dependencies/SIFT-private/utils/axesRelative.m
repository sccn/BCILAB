function [hAx hPan] = axesRelative(hParentAx, varargin)
% obtained from http://stackoverflow.com/questions/7113836/

    %# create panel exactly on top of parent axis
    s = warning('off', 'MATLAB:hg:ColorSpec_None');
    hPan = uipanel('Parent',get(hParentAx, 'Parent'), ...
        'BorderType','none', 'BackgroundColor','none', ...
        'Units',get(hParentAx,'Units'), 'Position',plotboxpos(hParentAx));
    warning(s)

    %# sync panel to always match parent axis position
    addlistener(handle(hParentAx), ...
        {'TightInset' 'Position' 'PlotBoxAspectRatio' 'DataAspectRatio'}, ...
        'PostSet',@(src,ev) set(hPan, 'Position',plotboxpos(hParentAx)) );

    %# create new axis under the newly created panel
    hAx = axes('Parent',hPan, varargin{:});
end