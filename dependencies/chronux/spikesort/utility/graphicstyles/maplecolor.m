function C = maplecolor(X,Y,Z)
%MAPLECOLOR        Generates truecolor data in the default Maple color scheme.
%   C = MAPLECOLOR(X,Y,Z) generates an M-by-N-by-3 matrix when each of X, Y,
%   and Z are M-by-N.  The resulting truecolor matrix maps [X,Y,Z] to
%   [R,G,B], scaled to soften the colors.

C = cat(3, X, Y, Z);

for slice = 1:3
    intense = C(:,:,slice);
  
    % Rescale each slice to [0.3, 0.9]
    intense = intense - min(intense(:));
    intense = intense ./ max(intense(:)) .* 0.6;
    intense = intense + 0.3;

    C(:,:,slice) = intense;
end
