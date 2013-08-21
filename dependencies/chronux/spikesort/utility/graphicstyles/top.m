function top
%TOP               Simplifies 3-D axes.
%   TOP deletes all lights in the current axes, turns off interpolated 
%   shading, and sets the view to overhead so that the surface looks like
%   an image.

delete(findobj(gca, 'Type', 'light'));
shading flat;    
lighting phong;     % otherwise you see seams in the surface.
view([0,90]);
