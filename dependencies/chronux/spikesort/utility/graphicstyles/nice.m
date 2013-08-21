function nice
%NICE              Custom plot prettification.
%   The function NICE is a macro for Matlab graphics settings, and its
%   behavior depends on the contents of the current axes.
%
%   If the plot contains only 2-D objects (i.e., objects with empty ZData)
%   NICE  sets the Y-limits depending on the data:
%               Data Range        Gives Limits
%                [ 0,1]           [-0.2 1.2]
%                [-1,1]           [-1.2 1.2]
%                [ 0,B]           [-0.2 1.2]*B
%                [-A,B]           [-1.2 1.2]*C, where C=MAX(ABS(A,B))
%   where the first condition that matches is used.
%   For data with range [P,Q], the X-limits are set to [P-Z, Q+Z], where
%   Z=0.02*MAX(ABS(P,Q)).
%
%   If the plot contains objects with non-empty ZData, 3-D visualization
%   aids are set, including: standard lighting conditions (existing lights
%   are removed), standard view angle for 3-D rotation and interpolated
%   shading (where possible).

childs = get(gca, 'Children');
if (length(childs) == 1),  childs = childs([1,1]);  end;  % ensure get(childs, prop) always returns cell

texts  = childs(strcmp(get(childs, 'Type'), 'text'));     % split into different types
lights = childs(strcmp(get(childs, 'Type'), 'light'));
images = childs(strcmp(get(childs, 'Type'), 'image'));
patchs = childs(strcmp(get(childs, 'Type'), 'patch'));
surfs  = childs(strcmp(get(childs, 'Type'), 'surface'));
lines  = childs(strcmp(get(childs, 'Type'), 'line'));
groups = childs(strcmp(get(childs, 'Type'), 'hggroup'));

surfpatch = [surfs(:); patchs(:)];
surfpatchline = [surfs(:); patchs(:); lines(:)];
surfpatchlinegroup = [surfs(:); patchs(:); lines(:); groups(:)];
%twoD = all(cellfun('isempty', get([lines,surfs,patchs], 'ZData')));
twoD = isequal(get(gca,'View'),[0 90]);

if (twoD)	
    allxdata = get(surfpatchlinegroup, 'XData');   allydata = get(surfpatchlinegroup, 'YData');   
	xdata = [];  ydata = [];
	for k = 1:length(surfpatchlinegroup)
		xdata = [xdata(:); allxdata{k}(:)];   ydata = [ydata(:); allydata{k}(:)];
	end
	
	ydatlims = minmax(ydata);
	maxabs    = max(abs([ydatlims,1]));
	if (ydatlims(1) < 0),   ylim = [-1.2, 1.2]*maxabs;
    else                    ylim = [-0.2, 1.2]*maxabs;
	end
	
	xdatlims = minmax(xdata);
	pad = max(0.02*abs(xdatlims));   pad = [-pad pad];
	xlim = xdatlims + pad;
	
	set(gca, 'YLim', ylim, 'XLim', xlim);
	zoom reset;

elseif (~isempty(surfpatch))
	
	axis tight;  set(gca,'CameraViewAngleMode', 'manual', 'CameraViewAngle', 9);

	interpable = ~cellfun('isempty', get(surfpatch, {'CData'}));
	set(surfpatch(interpable), 'FaceColor', 'interp', 'EdgeColor', 'none');

	material dull;	
	lighting gouraud;   % not as nice as phong lighting, but faster to render
	lightangle(90,60);  lightangle(270,10);
	
end