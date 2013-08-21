function h = mplot(t, matrix, varargin)
%MPLOT             Plot rows of a data matrix.
%   MPLOT(MATRIX) makes a line plot of the rows of an (M x N) matrix using
%   a single line object.  For large M, these plots are drawn more quickly
%   than PLOT(MATRIX') since only one object is created and that object
%   does not undergo front-to-back sorting.
%
%   MPLOT(X,MATRIX). where X is a length N vector, plots the rows of
%   MATRIX vs X.
% 
%   MPLOT(MATRIX, ...) passes additional arguments through to PLOT.
%
%   H = MPLOT(MATRIX) returns a handle to the line object.

%%%%% Parse arguments
if (nargin == 1),  matrix = t;  clear t;  end;
if ((nargin > 1) && ischar(matrix))  % MPLOT(MATRIX, ...) syntax
        varargin = {matrix; varargin{:}};
        matrix = t;   clear t;
end
if (~exist('t','var')),  t = 1:size(matrix,2);  end;

[M,N] = size(matrix);
[mm,nn] = size(t);
if (((mm == 1) && (nn == N)) || ((mm == N) && (nn == 1)))
    % do nothing here ...
else
    error('X must be a length N vector when MATRIX is M x N.');
end

%%%%% Create wraparound indices 
inds = repmat([t NaN]', [M,1]);
matrix = padmatrix(matrix, [0 1 0 0], NaN)';  % padding with NaNs keeps wraparounds from being drawn

%%%%% Do the plot.
h = plot(inds, matrix(:), varargin{:});

if (nargout < 1), clear h;  end;