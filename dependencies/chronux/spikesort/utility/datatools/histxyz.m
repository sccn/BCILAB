function [smoothdata,x_inds,y_inds,z_inds] = ...
	              histxyz(data_mtx, bins, smoothsd, levels, color, newfig)
%HISTXYZ           3-Dimensional Density Histogram.
%   HISTXYZ(DATA_MTX), where DATA_MTX is an (M x 3) matrix, visualizes 
%   the density of a 3-D scatter plot.  The rows generate a 3-D histogram
%   which is then smoothed with a 3-D Gaussian kernel.  This histogram is
%   treated as a scalar-valued volume and isocontours are found at a number
%   of density levels, each drawn as a transparent, colored surface.  The
%   coloring of the axes labels is (R,G,B) <=> (X,Y,Z).
%   Note that the fine details of the volume are simplified for graphics
%   efficiency beyond simple smoothing and will not necessarily be
%   accurate in details.
%
%   When called as [N,X,Y,Z] = HISTXYZ(...), the function returns the 3-D
%   histogram and associated indices rather than drawing the volume.
%
%   The surface is controlled by 4 optional arguments:
%       HISTXYZ(DATA_MTX, BINS, SMOOTHSD, LEVELS, COLOR, NEWFIG)
%
%   Any optional argument specified as [] will use its default value.
%   The definitions of the control arguments and there defaults are:
%      BINS            (default: 50) Number of bins used to discretize
%                       the data.
%      SMOOTHSD        (default: 1) The standard deviation in bin units
%                       for a Gaussian smoothing kernel.
%      LEVELS          (default: 10) A scalar value is taken as the #
%                       of evenly spaced density levels at which to
%                       draw a contour.  Vectors with values between
%                       [0..1] are treated as fractions of max density,
%                       e.g., [0.1 0.2 0.3] shows contours at 10%, 20%
%                       and 30% of max density.
%      COLOR           (default: evenly spaced intervals in the colormap
%                       are used for increasing density contours)  If
%                       specified as [R G B], it is taken to represent
%                       a single color at which all contours will be drawn.
%      NEWFIG          (default : 1) Set to 0 to draw on the current axes
%                       without changing camera properties.

% fraction of faces to be retained during patch simplification
patch_reduction = 0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (size(data_mtx,2) ~= 3),   error('Data must by [m x 3].');  end
if (nargin < 2 || isempty(bins)),      bins = 50;     end
if (nargin < 3 || isempty(smoothsd)),  smoothsd = 1;  end
if (nargin < 4 || isempty(levels)),    levels = 10;
else
    if (numel(levels) == 1)
        if (levels < 1),  error('If the ''levels argument is a scalar, it must be > 1.''');  end
    elseif ((ndims(levels) > 2) || (any(levels > 1) || any(levels < 0)))
        error('If the ''levels'' argument is a vector, it must be 2-D and it cannot have any values outside of [0..1].');
    end
end
if (nargin < 5 || isempty(color)),     color = colormap; end
if (nargin < 6 || isempty(newfig)),    newfig = 1;       end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rescale Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale the data to a range convenient for histogramming . . .
[data1,min1,max1] = rescale(data_mtx(:,1),1,bins);   data1 = round(data1);
[data2,min2,max2] = rescale(data_mtx(:,2),1,bins);   data2 = round(data2); 
[data3,min3,max3] = rescale(data_mtx(:,3),1,bins);   data3 = round(data3);

% Construct index vectors from the original ranges of the data.
x_inds = linspace(min1,max1,bins);
y_inds = linspace(min2,max2,bins);
z_inds = linspace(min3,max3,bins);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the 3-D scatter density
% 'sparse' only works in 2-D so we combine the first 2 dimensions ...
counts = full(sparse(data2+(bins)*(data1-1), data3, 1, bins^2, bins));
counts = reshape(counts,[bins,bins,bins]);  % . . . and now separate them.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout > 0),  return;  end;
% First we smooth:
kernelsize = ceil(smoothsd*8 + 1);
smoothdata = smooth3f(counts, 'gaussian', kernelsize, smoothsd);

% Now determine which contour levels we'll want to draw.
maxdens = max(smoothdata(:));
onedens = max(max(max(gausskernel(kernelsize*ones(1,3), smoothsd))));
if (numel(levels) == 1)
	showlevels = linspace(onedens,maxdens,levels+1);
	showlevels = showlevels(1:end-1);
else
	showlevels = round(sort(levels) * maxdens);
end

% Set up evenly spaced colors.
cmap = color;
cind = floor(linspace(1,size(cmap,1),length(showlevels)));

% Prepare figure and prettify.
if (newfig)
	figure;
	grid on;   set(gca,'box','on');   daspect([1,1,1]);   view(3);
end
axhand = gca;

% Next we actually draw the contours.
for contour = 1:length(showlevels)
	fv = isosurface(x_inds,y_inds,z_inds,smoothdata,showlevels(contour));
	fv = reducepatch(fv, patch_reduction);
	p = patch(fv);
	
	set(p, 'FaceColor', cmap(cind(contour),:), 'EdgeColor', 'none');
	% The alpha decreases inversely (to countour #) to emulate increasing density.
	set(p, 'AlphaDataMapping', 'none', 'FaceAlpha', 1/(length(showlevels)));
	
	drawnow;
end

% Prepare the axes for rotating and add some lighting 
if (newfig)
	axis tight;
	set(axhand,'CameraViewAngleMode','manual', 'CameraViewAngle', 10);
end
material dull;  
lighting gouraud;   lightangle(90,60); 	lightangle(270,10);

% Finally, color the tick marks . . . 
set(axhand, 'Xcolor', [0.6 0 0]);
set(axhand, 'Ycolor', [0 0.6 0]);
set(axhand, 'Zcolor', [0 0 0.6]);

% . . . and make sure that nothing is returned.
clear smoothdata x_inds y_inds z_inds
