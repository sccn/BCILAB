function map = colormap_visualization(map, c)
%colormap_visualization: produce visualisations for evaluating colormaps.
% Example:
%  colormap_vis(jet(9))
%
% See also: colormap_investigation, colormap_optimization, bipolar
% and real2rgb, available from MATLAB Central Fileexchange ID 23342:
%  http://www.mathworks.com/matlabcentral/fileexchange/23342

% Copyright 2009 Ged Ridgway at gmail com

if nargin < 2
    c = 3; % everything
end

%% Raw colormap
if c == 1 || c == 3
    figure
    membrane(3, 25); axis equal; colormap(map); colorbar
    drawnow

    %% Interpolated to 256 levels with Oliver Woodford's real2rgb (#23342)
    tmp = membrane(3, 25);
    [tmp tmp map] = real2rgb(tmp, map); %#ok tmp unused
    figure
    membrane(3, 25); axis equal; colormap(map); colorbar
    drawnow
end

%% Converted to grayscale (luminance) using Gimp's formula
% From http://gimp-savvy.com/BOOK/index.html?node54.html
% See colormap_optimization for more details
if c == 2 || c == 3
    coef = [0.3 0.59 0.11]';
    map = map*[coef coef coef];
    figure
    subplot(1,2,1)
    membrane(3, 25); axis equal; colormap(map); colorbar

    subplot(1,2,2)
    plot(map*coef);
    xlabel('Index into colormap');
    ylabel('Luminance value');
    ylim([0 1])
    drawnow
end

if nargout == 0, clear map, end
