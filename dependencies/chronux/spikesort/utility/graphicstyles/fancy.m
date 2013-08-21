function fancy(select)
%FANCY             Collection of graphics test figures.
%   fancy(1)    Sphere with a cloudy (transparency-mapped) ring
%   fancy(2)    World map, mapped onto a sphere, with transparent oceans
%   fancy(3)    3-D relief model of a US penny, with colored lighting
%
% Each of the figures will briefly animate after creation.

switch(select),
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case(1),  %%% sphere w/ ring
figure
resol  = 50;    % # pts for sphere & theta pts for ring
resol2 = 25;    % radial pts for ring
blur   = 13;    % feature blurring on sphere
blur2  = 5;     % alpha blurring on ring

set(gcf,'DoubleBuffer','on', 'Color',[0 0 0]);
colormap(autumn(256));

[X1,Y1,Z1] = sphere(resol-1);                   % make inner sphere
C1 = conv2(Z1(reshape(randperm(resol.^2), resol, resol)), gausswin(blur) * gausswin(blur)', 'same');
C1 = staggermatrix(C1);      % creates a 'swept' look by breaking latitudinal blurring only

[TH, R] = meshgrid(linspace(0,2*pi,resol), linspace(1.5, 2.5, resol2));  % make ring
[X2,Y2] = pol2cart(TH', R');
Z2 = zeros(resol, resol2);
A2 = conv2(rand(resol,resol2), gausswin(blur2) * gausswin(blur2)','same');  % random alpha map, blurred
A2 = staggermatrix(A2')';    % good for rings too, since it breaks the radial blurring
A2(1,:) = A2(end,:);         % prevents a seam by allowing interpolation from both sides  

h1 = surf(X1,Y1,Z1,C1);          % draw inner sphere
set(h1, 'Ambient', 0.6);

hold on                                         % draw ring
h2 = surf(X2,Y2,Z2);
set(h2, 'FaceAlpha','interp', 'AlphaData', A2);

shading interp; lighting gouraud; material dull;
axis equal; axis off;
set(gca,'CameraViewAngleMode','manual', 'CameraViewAngle', 7);
lightangle(90,40);

set(h2, 'FaceColor', [0.7, 0.7, 0.7]);  % ring should be plain gray in color

% Animate
view([0,15]);
for steps = 1:72, 
	camorbit(10,0); drawnow;
	pause(0.010); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case(2),  %%% display globe on a sphere with transparent oceans
figure;  load topo;
topo = padmatrix(topo, [0 0 0 10], -3500);  % the topo map leaves off the top of the globe ...

[X,Y,Z] = sphere(180);
h1 = surf(X,Y,Z);

set(h1, 'FaceColor', 'texture', 'EdgeColor', 'none', 'FaceAlpha', 'texture', 'AlphaDataMapping', 'none');
set(h1, 'CData', topo, 'AlphaData', (double(topo > 0) .* 0.6) + 0.4);

colormap(topomap1);
axis equal;
set(gca, 'CameraViewAngleMode', 'manual');
axis off;

lt1 = lightangle(0,0);  set(lt1, 'Color', [1, 1, 1]);

% Animate
view([0,5]);
for step = 1:24, 
	camorbit(15, 0); drawnow; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case(3),  %%% create a penny surface model with lights, etc
figure;  load penny;
P = fliplr(P);
coords = -63:64;  [X,Y] = meshgrid(coords);
coppercolor = copper(2);  coppercolor = coppercolor(2,:);

mask = double(sqrt((X-0).^2 + (Y-12).^2) < 80);
mask(~mask) = NaN;
surf(X,Y,mask.*P, 'FaceColor', coppercolor, 'EdgeColor', 'none');

axis equal;  axis off;
set(gca,'CameraViewAngleMode','manual', 'CameraViewAngle', 9, 'DataAspectRatio', [1 1 30]);

material metal; lighting gouraud;
lt1 = lightangle(90, 10);    set(lt1, 'Color', 'r');
lt2 = lightangle(180, 30);   set(lt2, 'Color', 'g');
lt3 = lightangle(270, 60);   set(lt3, 'Color', 'b');
lt4 = lightangle(0, -90);    set(lt4, 'Color', 'w');

% Animate
view([160,50]);
for step = 1:71
	camorbit(3, 20); drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otherwise,
	error('Unknown selection.');
end

%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%

function A = staggermatrix(A)
% STAGGERMATRIX  Randomly circularly shifts the rows of an input matrix.
%    B = STAGGERMATRIX(A) circularly shifts each row of the matrix A by a
%    random (integer) amount.

[M, N] = size(A);
for row = 1:M
    shifts = randperm(N);
    A(row, :) = circshift(A(row, :), [1, shifts(1)]);
end
