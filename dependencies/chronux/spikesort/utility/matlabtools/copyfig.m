function target = copyfig(source)
%COPYFIG           Duplicate all children of a figure.
%   COPYFIG(FIG) makes a copy of the figure specified by the handle FIG in
%   a new figure.  The copied figure includes all children of FIG.
%
%   FIG2 = COPYAXES(FIG) returns a handle to the copied FIGURE.
%
%   See also COPYOBJ, COPYAXES.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~ishandle(source) || ~strcmp(get(source,'Type'), 'figure'))
	error('Input must be the handle to an existing figure.');
end
children = get(source, 'Children')';
possource = get(source, 'Position');
clrmap = get(source, 'Colormap');
renderer = get(source, 'Renderer');

%%%%%%%%%%%%%%%%%%%%%%%%% Create Target Figure %%%%%%%%%%%%%%%%%%%%%%%%%
target = figure;
set(target, 'Position', [50 50 possource(3:4)]);  % copy orig size but not location
set(target, 'Colormap', clrmap, 'Renderer', renderer);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the Copy %%%%%%%%%%%%%%%%%%%%%%%%%%%%
copyobj(children, target);
if (nargout == 0),  clear target;  end;