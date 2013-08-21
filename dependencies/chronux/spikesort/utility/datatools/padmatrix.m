function padded = padmatrix(input, wherepad, padval)
%PADMATRIX         Pad a matrix.
%   PADDED = PADMATRIX(INPUT, WHEREPAD) adds zero rows/columns to the
%   matrix INPUT.  WHEREPAD is a four element vector specifying the number
%   of pad rows/cols beyond [left col, right col, top row, bottom row].
%
%   PADDED = PADMATRIX(INPUT, WHEREPAD, PADVAL) pads with PADVAL elements
%   rather than zeros.
%
%   Simplified version of Matlab's Image Processing Toolbox' PADARRAY.

%%%%%%%%%% ARGUMENT CHECKING
if (nargin < 3),  padval = 0;  end;

if (ndims(input) > 2)
    error('PADMATRIX currently does not handle arrays of dimension greater than 2.');
elseif ((numel(wherepad) ~= 4) && (length(wherepad) ~= 4))
	error('Second argument must be a four-element vector.');
elseif (numel(padval) > 1)
	error('Third argument, when specified, must be a scalar.');
end

[M,N] = size(input);
newN = N + sum(wherepad(1:2));

%%%%%%%%%% DO THE PAD
% columns first
padleft = repmat(padval, [M, wherepad(1)]);
padright = repmat(padval, [M, wherepad(2)]);
padded = [padleft input padright];

% then rows
padtop = repmat(padval, [wherepad(3), newN]);
padbottom = repmat(padval, [wherepad(4), newN]);
padded = [padtop; padded; padbottom];

