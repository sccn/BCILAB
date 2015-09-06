function new = colormap_optimization(cmap)
%colormap_optimisation: minimally modify colormap to satisfy constraints.
% In this case, the constraints are for linear grayscale luminance and a
% neutral mid-gray central color.
%
% Example:
%  sprung = colormap_optimization(spring(5));
%  colormap_visualization(spring(5))
%  colormap_visualization(sprung)
%
% See also: colormap_investigation, colormap_visualization, bipolar

% Copyright 2009 Ged Ridgway at gmail com

%% Conversion from RGB to grayscale brightness/luminance

% From http://www.poynton.com/ColorFAQ.html (Q.9)
% coef = [0.2126 0.7152 0.0722]'

% From rgb2gray:
% T = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);
% coef = T(1,:)';

% From  http://gimp-savvy.com/BOOK/index.html?node54.html
coef = [0.3 0.59 0.11]';
% Note, the Gimp coef is almost identical to that from rgb2gray, and these
% seem nicer in a print('-dps') grayscale copy than the Poynton version.

%% Constraints
m = size(cmap, 1);
mid = zeros(1, m);
mid(round((m+1)/2)) = 1;
Aeq = [
    kron(coef', eye(m)) % linearity in luminance
    kron(eye(2,3), mid) % centrality
    ];
l = 0.1; % lower luminance (avoid black, to give contrast with text/lines)
u = 0.9; % upper luminance (avoid white, to give contrast with background)
beq = [
    linspace(l, u, m)'  % linearity in luminance
    ones(2, 1) / 2      % centrality
    ];

%% Optimisation (using MATLAB optimisation toolbox)
orig = warning('off', 'optim:lsqlin:LinConstraints');
x = lsqlin(eye(3*m), cmap(:), [], [], Aeq, beq, zeros(3*m, 1), ones(3*m, 1));
x = max(min(x, 1), 0); % enforce bounds, since lsqlin can give e.g. -1e-18
new = reshape(x, m, 3);
warning(orig);
