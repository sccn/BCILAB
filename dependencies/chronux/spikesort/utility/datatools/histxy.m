function [counts,x_inds,y_inds] = histxy(x, y, varargin)
%HISTXY            2-Dimensional Density Histogram.
%   COUNTS = HISTXY(X,Y,D), where X and Y are matrices of the same size,
%   returns a D x D matrix COUNTS containing the numbers of points in a
%   binned scatter plot of the entries of X vs the corresponding entries
%   of Y.  The elements of X and Y are each independently discretized
%   into D evenly-spaced bin centers, offset such that the extrema are
%   are bin centers.  If D is not specified (or is empty) it defaults
%   to 100.
%
%   COUNTS = HISTXY(X,Y,[D1,D2]) produces a D2 x D1 matrix COUNTS, with
%   the corresponding binning of X and Y.
%
%   [COUNTS,X_INDS,Y_INDS] = HISTXY(X,Y,D) also returns the bin centers.
%   The density is then visualized with IMAGESC(X_INDS,Y_INDS,COUNTS).
%
%   [...] = HISTXT(X,Y,D,BANDWIDTH), for scalar BANDWIDTH, convolves
%   COUNTS with an isotropic 2-D Gaussian kernel of standard deviation
%   BANDWIDTH bins.  If BANDWIDTH is a two element vector, an anisotropic
%   Gaussian is used with column standard deviation BANDWIDTH(1) bins and
%   row standard deviation BANDWIDTH(2) bins.
%
%   [...] = HISTXY(...,'log') uses the log of the counts (0's yield -Inf).
%
%   HISTXY(...) without output arguments produces an image of the density.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = [100,100];  logflag = 0;   blur = 0;      % defaults

x = x(:);  M = length(x);
y = y(:);  N = length(y);
if ((M~=N) || ~isnumeric(x) || ~isnumeric(y) || any(isinf(x(:))) || any(isinf(y(:))))
	error('First 2 input arguments must be numeric matrices of equal size with no +/- Inf elements.');
end
mask = (isnan(x) | isnan(y));   nanflag = any(mask(:));
if (nanflag),  warning('HISTXY:NaNflag', 'NaN elements will be ignored.');  end

if (length(varargin) > 0)
	tail = varargin{end};
	if (ischar(tail) && strcmpi(tail,'log'))      % If the last arg was 'log' ...
		varargin = varargin(1:(end-1));           % ... chomp it and set a flag.
		logflag = 1;
	end
	if (length(varargin) > 1)                     % If two args left, ...
		tail = varargin{end};                     % ... try to chomp bandwidth.
		if (isnumeric(tail) && length(tail)<=2)
			varargin = varargin(1:(end-1));
			if (length(tail) == 2),  sigma = tail;
			elseif (length(tail) == 1), sigma = [tail, tail];
            else   error('Bandwidth can not be empty.');
			end
			blur = 1;
			sigma(sigma == 0) = 0.01;  % equiv to 0 b/c 100 stds till next bin => below realmin
		end
	end
	if (length(varargin) > 0)
		tail = varargin{end};
		if (isnumeric(tail) && length(tail)<=2)   % If remaining arg was a len<=2 vector, ...
			varargin = varargin(1:end-1);         % ... chomp it and set the bin count.
			if (length(tail) == 2),  D = tail;
			elseif (length(tail) == 1),  D = [tail, tail];
			end  % if tail was empty, use default ...
		end
	end
	if (length(varargin) > 0),  error('Unknown syntax.');  end;
end

if (~isequal(D,round(D))), error('Number of bins D must be integer-valued.');  end;
if (any(D==1)),  error('Number of bins D must be greater than 1.');  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rescale Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separately scale the x and y data
[x,min1,max1] = rescale(x,1,D(1));   x = round(x);
[y,min2,max2] = rescale(y,1,D(2));   y = round(y);

% Bin centers are equally spaced over the range of the unscaled data
x_inds = linspace(min1,max1,D(1));
y_inds = linspace(min2,max2,D(2));

% Mask NaNs
if (nanflag),  D=D+1;  x(mask)=D(1);  y(mask)=D(2);  end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counts = CORE_histxy(x,y,D(1),D(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clean Up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nanflag),  counts(end,:) = []; counts(:,end) = [];   D = D-1;  end;

if (blur)   % use separable kernels for convolution with least odd length  >= num bins
	x_kernel = exp(-((-floor(D(1)/2):floor(D(1)/2)).^2)./(2*sigma(1).^2));
	x_kernel = x_kernel ./ sum(x_kernel);
	y_kernel = exp(-((-floor(D(2)/2):floor(D(2)/2)).^2)./(2*sigma(2).^2))';
	y_kernel = y_kernel ./ sum(y_kernel);
	
	counts = conv2(y_kernel, x_kernel, counts, 'same');
end

if (logflag), o=warning('MATLAB:log:logOfZero', 'off');  counts=log(counts);  warning(o);  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout > 0),  return;  end;
imagesc(x_inds, y_inds, counts);   axis xy;
clear counts x_inds y_inds  % clear these so nothing is dumped to output
