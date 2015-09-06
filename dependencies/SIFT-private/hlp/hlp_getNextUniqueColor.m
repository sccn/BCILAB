function color = hlp_getNextUniqueColor(varargin)
% Return the next available color from a color map. This function will
% cycle circularly through the colormap
%
% Inputs:
%   colorMap: A [num_entries x 3] colormap matrix. This need only be
%             supplied once, upon which it will be stored persistently and
%             automatically indexed until a new colormap is provided.
%   'reset':  if supplied as the final argument, the index will be reset
%
% Example:
%
% jet(4)
%
%      0     0     1
%      0     1     1
%      1     1     0
%      1     0     0
%
% color = hlp_getNextUniqueColor(jet(4))   % get the next color
%
%      0     0     1
%
% color = hlp_getNextUniqueColor            % get the next color
%
%      0     1     1
%
% color = hlp_getNextUniqueColor('reset')   % reset the index
%
%      0     0     1
%
% color = hlp_getNextUniqueColor(jet(7))    % change the colormap (resets the index)
%
%      0     0     1

persistent index;
persistent cmap;

doReset = false;
colorMap = cmap;

if (nargin==0 || (nargin==1 && any(strcmpi(varargin,'reset')))) ...
    && isempty(colorMap)
    % initialize with a default colormap
    colorMap = hsv(10);
end

if nargin>0
    if ~ischar(varargin{1})
        % user supplied a colormap
        colorMap = varargin{1};
    end
    
    if ~isequal(cmap,colorMap)
        % ..colormap is new, update
        doReset = true;
        cmap = colorMap;
    end
end

doReset = doReset | any(strcmpi(varargin,'reset'));

% (re-)initialize index
if doReset || isempty(index)
    index = 0;
end


index = index + 1;

% ring buffer
index = mod(index-1,size(colorMap,1))+1;

color = colorMap(index,:);