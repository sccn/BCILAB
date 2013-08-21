function raster = rasterize(events, Fs, T)
%RASTERIZE         Converts a list of events into a binary time series.
%   RASTER = RASTERIZE(EVENTS) takes a vector EVENTS containing a list of
%   indices and returns a binary vector RASTER with ones in locations
%   given by EVENTS.  E.g., RASTERIZE([1 5]) => [1 0 0 0 1].
% 
%   RASTER = RASTERIZE(EVENTS,FS) specifies a sampling frequency FS;
%   values in EVENTS are multiplied by FS and rounded to obtain indices
%   into RASTER; e.g., RASTERIZE([0.1, 0.47], 10) => [1 0 0 0 1].
%   RASTERIZE(EVENTS) is equivalent to RASTERIZE(EVENTS,1).
%
%   RASTER = RASTERIZE(EVENTS,FS,T) specifies a total time T, so that
%   RASTER is length FS*T (events outside of this interval are ignored).
%   The default value is the maximum of EVENTS.
%
%   The array returned for RASTER is a row vector of type LOGICAL.

%%%%% Deal with Fs & T ...
if (nargin > 1),  events = round(events*Fs);
else              events = round(events);
end

if (nargin > 2),  len = Fs*prod(T);  
else              len = max(events);  
end

numevents = length(events);
events(events>len) = [];
if (length(events) < numevents), warning('Utility:raster_overflow', 'Some events exceed raster grid.'); end;

numevents = length(events);
events(events<=0)  = [];
if (length(events) < numevents), warning('Utility:raster_overflow', 'Some events preceed raster grid.'); end;

%%%%% Make the logical raster
raster = repmat(false, 1, round(len));
raster(events) = 1;
