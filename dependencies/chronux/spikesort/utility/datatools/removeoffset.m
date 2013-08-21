function Y = removeoffset(X, dim, mode, dc)
%REMOVEOFFSET      Remove DC offset using various techniques.
%   Y = REMOVEOFFSET(X) subtracts the mean value from the data in vector X
%   and returns it in vector Y.  If X is a matrix, REMOVEOFFSET removes
%   the mean value from each column.
%
%   Y = REMOVEOFFSET(X, DIM) removes the mean along dimension DIM for the
%   N-D matrix X.  DIM can be the empty matrix [], in which case it
%   defaults to 2 for row vectors and 1 for all other arrays.
%
%   Y = REMOVEOFFSET(X, DIM, 'median') subtracts the median.
%   Y = REMOVEOFFSET(X, DIM, 'mean')   subtracts the mean (default).
%   Y = REMOVEOFFSET(X, DIM, 'local')  subtracts a local 3x3 average.
%
%   Y = REMOVEOFFSET(X, DIM, MODE, DC) uses the matrix DC to compute the
%   offset, which is then removed from X to give Y.  If MODE is 'mean or
%   'median', DC must have the same size as X in all dimensions except for
%   DIM.  If MODE is 'local', DC must be the same size as X in all
%   dimensions.
%   
%   See also DETREND.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 4),  dc = X;  end;
if ((nargin < 3) || (isempty(mode))),  mode = 'mean';  end
if ((nargin < 2) || (isempty(dim)))
	if (isvectord(X)>1),  dim = 2;  else  dim = 1;  end;
end
if (ischar(dim)),  error('Second argument must be numeric.');  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute Offset %%%%%%%%%%%%%%%%%%%%%%%%%%%
% need this to match the computed offsets dimensions to X
sz = size(X);  rep = ones(size(sz));  rep(dim) = sz(dim);

switch (mode),
	case 'mean',    offset = repmat(mean(dc, dim), rep);
	case 'median',  offset = repmat(median(dc, dim), rep);
	case 'local',   offset = conv2(dc, ones(3)/9, 'same');
	otherwise,      error('Invalid mode.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Remove Offset %%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = X - offset;