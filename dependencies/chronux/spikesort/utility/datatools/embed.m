function embedded = embed(timeseries, lag, D)
%EMBED             Time-series lag embedding.
%   EMBEDDED = EMBED(TIMESERIES, LAG_SAMPLES, DIMENSION) takes a scalar
%   (M x 1) time series and returns a vector time series of dimension
%   ((M - lag_samples*(dimension-1)) x dimension) time series by lag
%   embedding the time series.  E.g.,
%             embed([x(1); x(2); x(3); ... ; x(N)], 1, 3)
%   returns
%             [x(1), x(2), x(3); x(2), x(3), x(4); ... ; x(N-1), x(N)]
%
%   NOTE: A matrix input 'timeseries' is treated as 'timeseries(:)'

if (nargin < 3)
	error('Incorrect number of arguments.');
end

% Column-ize the time series.
timeseries = timeseries(:);

M = length(timeseries);  % length of input time series
N = M - lag * (D-1);     % length of output time series

% Now do the embedding efficiently by indexing into the original timeseries
embedded = timeseries(repmat((1:N)', 1, D) + repmat(lag * (0:(D-1)), N, 1));

