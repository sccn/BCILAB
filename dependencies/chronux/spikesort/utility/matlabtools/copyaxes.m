function target = copyaxes(source)
%COPYAXES          Duplicate axes in a new figure.
%   COPYAXES(AX) makes a copy of the axes specified by the handle AX in a
%   new figure.  The copied axes are identical to the original except that
%   their position is modified to fill the new parent figure.  This can be
%   used to enlarge an axes in a figure subplot.
%
%   AX2 = COPYAXES(AX) returns a handle to the copied axes.
%
%   See also COPYOBJ, COPYFIG.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~ishandle(source) || ~strcmp(get(source,'Type'), 'axes'))
	error('Input must be the handle to an existing axes.');
end

%%%%%%%%%%%%%%%%%%%%%%%%% Create Target Axes %%%%%%%%%%%%%%%%%%%%%%%%%
hfig = figure;
position = get(gca, 'Position');  % default axes position for single subplot
delete(gca);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the Copy %%%%%%%%%%%%%%%%%%%%%%%%%%%%
target = copyobj(source, hfig);
set(target, 'Position', position);

if (nargout == 0),  clear target;  end;