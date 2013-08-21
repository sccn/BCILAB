function meandist = meandist_estim(data, use_data)

% MEANDIST_ESTIM  Estimate the average distance between vectors.
%    MEANDIST = MEANDIST_ESTIM(DATA) takes an (M x N) matrix DATA and returns
%    an estimate of the average Euclidean distance between the row vectors
%    of DATA.  The estimate is made by averaging the pairwise distances for
%    500 randomly selected rows.
% 
%    MEANDIST = MEANDIST_ESTIM(DATA, USE_DATA), when USE_DATA is a scalar, uses
%    min(M, USE_DATA) rows for the estimate.  If USE_DATA is a vector, it is
%    taken as a set of indices into the rows of DATA and the distance estimate
%    is made with the specified rows.

%%%%%%%%%% CONSTANTS
default_numrows = 800;
M = size(data, 1);

%%%%%%%%%% DEFAULTS & ARGUMENT CHECKING 
if ((nargin == 2) && (numel(use_data) > 1))
	try
		data(use_data,:);    % lazy way of checking the indices
	catch
		error('The USE_DATA vector contains invalid indices.');
	end
else
	if (nargin == 1)
		use_data = min(M, default_numrows);
	elseif (use_data > M)
		error('USE_DATA can not contain a scalar higher than the number of rows in DATA.');
	end
	inds = randperm(M);   % randomly shuffle row indices ...
	use_data = inds(1:use_data);      % ... and take the number of indices we're going to use
end

%%%%%%%%%% MEAN DISTANCE ESTIMATION
dists = pairdist(data(use_data,:));
meandist = mean(dists(triu(logical(dists ~= 0))));

