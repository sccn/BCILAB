function h = topoplot_grid(W,chanlocs,varargin)
% like topoplot(), but takes an array of topographic weights and plots them in a 2d grid
% topoplot_grid(Weights,Chanlocs,Options...)
%
% In:
%   Weights    : [#channels x #maps] matrix of topographic maps (note: if one of the dimenions matches
%                the number of supplied chanlocs and the other doesn't, the orientation of the matrix
%                will be fixed if wrong)
%
%   Chanlocs   : EEGLAB chanlocs structure or location file)
%
%   Options... : optional name-value pairs, as in topoplot(); additional options are:
%                'titles' : cell array of titles; one per map

% process optional inputs
k = find(strcmp(varargin(1:2:end),'titles'),1,'last');
if isempty(k)
    titles = {};
else
    titles = varargin{2*k};
    varargin = varargin([1:2*k-2 2*k+1:end]);
end

k = find(strcmp(varargin(1:2:end),'rows'),1,'last');
if ~isempty(k)
    rows = varargin{2*k};
    varargin = varargin([1:2*k-2 2*k+1:end]);
end
k = find(strcmp(varargin(1:2:end),'cols'),1,'last');
if ~isempty(k)
    cols = varargin{2*k};
    varargin = varargin([1:2*k-2 2*k+1:end]);
end
k = find(strcmp(varargin(1:2:end),'scales'),1,'last');
if ~isempty(k)
    scales = varargin{2*k};
    varargin = varargin([1:2*k-2 2*k+1:end]);
end
k = find(strcmp(varargin(1:2:end),'handles'),1,'last');
if ~isempty(k)
    handles = varargin{2*k};
    varargin = varargin([1:2*k-2 2*k+1:end]);
end

if isnumeric(titles)
    titles = cellfun(@num2str,num2cell(titles),'UniformOutput',false); end

if ischar(chanlocs)
    chanlocs = readlocs(chanlocs); end

% flip orientation if necessary
if size(W,1) ~= length(chanlocs) && size(W,2) == length(chanlocs)
    W = W'; end

maps = size(W,2);
if ~exist('cols','var')
    cols = ceil(sqrt(maps)); end
if ~exist('rows','var')
    rows = ceil(maps/cols); end

if isstruct(chanlocs) && ~isfield(chanlocs(1),'sph_theta')
    chanlocs = convertlocs(chanlocs,'cart2all');
end

% for each map...
for p=1:maps
    if exist('handles','var')
       [dummy grid] = topoplot(W(:,p),chanlocs,varargin{:},'noplot','on');
       set(handles{p},'Cdata',grid);
    else
        subplot(cols,rows,p);  
        h{p} = topoplot(W(:,p),chanlocs,varargin{:});
    end
    if exist('scales','var')
        camzoom(scales(p)); end
    if ~isempty(titles)
        try,title(titles{p},'Interpreter','none'),catch,end; end
    if mod(p,rows)==0
        drawnow; end
end
