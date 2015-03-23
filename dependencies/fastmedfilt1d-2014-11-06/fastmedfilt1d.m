function m = fastmedfilt1d(x, W, xic, xfc)
% Very fast implementation of the classical 1D running median filter of
% window size W (odd). Uses a quick selection method to obtain the median
% of each sliding window. Initial and final condition vectors can also be
% specified. Avoids the need for a large sort buffer so handles very long
% window sizes efficiently.
%
% Usage:
% m = fastmedfilt1d(x, W, xic, xfc)
%
% Input arguments:
% - x          Input signal
% - W          Window size (W = 3 if not specified)
% - xic        Initial condition (all zeros if not specified)
% - xfc        Final condition (all zeros if not specified)
%
% Output arguments:
% - m          Median filtered output signal
%
% (c) Max Little, 2010. If you use this code, please cite:
% Little, M.A. and Jones, N.S. (2010),
% "Sparse Bayesian Step-Filtering for High-Throughput Analysis of Molecular
% Machine Dynamics"
% in Proceedings of ICASSP 2010, IEEE Publishers: Dallas, USA.

% Input parameter case handling
x = x(:);
N = length(x);
error(nargchk(1,4,nargin));
if (nargin == 1)
    W = 3;
    xic = zeros(1,1);
    xfc = xic;
end

% Ensure that W is odd
W2 = floor(W/2);
W  = 2*W2 + 1;

if (nargin == 3)
    xfc = zeros(W2,1);
end
if (nargin == 2)
    xic = zeros(W2,1);
    xfc = xic;
end

% Input parameter checking
if (W < 1)
    error('Window size W must be greater than zero');
elseif (W >= floor(N/2))
    error('Window size W must be less than half the length of x');
end
if ((length(xfc) ~= W2) || (length(xic) ~= W2))
    error('Initial/final conditions must be a vector of length floor(W/2)');
end

xic = xic(:);
xfc = xfc(:);

% Call fast MEX core routine
m = fastmedfilt1d_core(x, xic, xfc, W2);
