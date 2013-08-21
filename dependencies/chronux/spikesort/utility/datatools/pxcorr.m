function [C, lags] = pxcorr(x, varargin)
%PXCORR            Efficient cross-correlation for point process data.
%   C = PXCORR(A,B,Fs) returns the cross-correlation between A and B at
%   sampling rate Fs for the point process whose times are given in the
%   length N(>1) vector A and the length M (>1) vector B.  C will be a row
%   vector with length given by
%                2*Fs*max(max(A)-min(B),max(B)-min(A))+1
%   i.e., twice the maximum time difference between events in A and B
%   measured in units of 1/Fs (plus 1 for 0 lag).  The vector C estimates
%                       C(t) = E[A(x-t)*B(x)]
%   The event times listed in A and B are assumed to be sorted and will be
%   rounded to bins of duration 1/Fs. 
%
%   The computation here is O(N*J), where 0<=J<=M is the mean number of
%   events in the vector B that are within MAXLAG (see below) of events in
%   the vector A (when MAXLAG is not specified, J=M).  Matlab's XCORR is
%   O(K log K), where K = max(max(A),max(B))*Fs -- note that if MAXLAG is
%   large and the mean interevent interval is small, XCORR can be faster.
%
%   C = PXCORR(A,Fs) returns the auto-correlation of the events listed in A 
%   and is equivalent to PXCORR(A,A,Fs)
%
%   C = PXCORR(...,MAXLAG) only considers lags upto MAXLAG; the length of C
%   will then be 2*Fs*MAXLAG + 1.
%
%   PXCORR(...,'sort') indicates the inputs are not yet sorted.
%
%   [C,LAGS] = PXCORR(...) also returns a vector of lag indices.


%%%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<2), error('Not enough arguments.');  end;
if (length(x)<=1), error('Point processes must consist of >1 event.'); end;
if (~isvectord(x)), error('Matrix data are not supported.');  end;

presorted = 1;   % already sorted events?
if (strcmpi(varargin{end},'sort')),
    presorted = 0;
    varargin = varargin(1:end-1);   
	if(isempty(varargin)), error('Not enough arguments.'); end;
end;

if (length(varargin{1}) == 1), y = x;    % make auto-corr case look like X
else
    y = varargin{1};
    varargin = varargin(2:end);
	if(isempty(varargin)), error('Not enough arguments.'); end;
end
if (length(y)<=1),  error('Point processes must consist of >1 event.');  end;
if (~isvectord(x)), error('Matrix data are not supported.');  end;

x = x(:)';  y = y(:)';    % enforce row vectors
if (~presorted),  x = sort(x);  y = sort(y);  end;   % only do this if needed

nargs = length(varargin);
Fs = varargin{1};           
if(length(Fs)~=1), error('Fs must be a scalar.'); end;

if (nargs == 1)   % limiting the lag saves time
    maxlag = max(x(end)-y(1),y(end)-x(1));
elseif (nargs == 2)
    maxlag = varargin{2};           
	if(length(maxlag)~=1 || isinf(maxlag)), error('MAXLAG must be a finite scalar.'); end;
elseif ((nargs == 3) && (length(varargin{3})<=1))
    error('Point processes must consist of > 1 event.');
else
    error('Invalid syntax.');
end

x = round(x*Fs);  y = round(y*Fs);
maxlag = ceil(maxlag * Fs);

if (nargout == 2),  lags = (-maxlag:maxlag)./Fs;  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correlate %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = CORE_pxcorr(x,y,maxlag);