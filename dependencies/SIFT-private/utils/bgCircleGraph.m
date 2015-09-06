function hbg = bgCircleGraph(bg,NodeSizes,varargin)
% (optional) NodeSizes should be an array specifying the radius of each node
% as a fraction (0 to 1) of the radius of the circlegraph

if nargin<2 
    NodeSizes = [];
elseif numel(NodeSizes)==1
    NodeSizes = NodeSizes*ones(1,length(bg.Nodes));
end

% I create a random graph with 28 nodes:
% bg = biograph(28);

% Display it in the GUI and get a handle 
% back to work with it (this will use by default 
% the hierarchical layout):
hbg = view(bg);
% Get an idea of what are the actual extents of 
% the page in the GUI by looking at all the current 
% node positions:
page_Size = max(cell2mat(arrayfun(@(x) get(x,'Position'),...
            get(hbg,'Nodes'),'Uniform',false)));
% Place the nodes in a circular layout, in my graph 
% page_size was around [1000,1000], so I will select a 
% center at [500 500] and a radius of 300, this will keep 
% the layout approximately in the same scale and I will 
% not have to manually change the size of the nodes or 
% the fonts:
radius = 0.6*(min(page_Size)/2);
center = [page_Size(1)/2 page_Size(2)/2];
numNodes = length(hbg.Nodes);
for i = 1:numNodes
set(hbg.Nodes(i),'Position',...
   [center(1)+radius.*sin((i*2*pi/numNodes)),...
    center(2)+radius.*cos((i*2*pi/numNodes))],'shape','circle')

    if ~isempty(NodeSizes)
        sz = 2*radius*NodeSizes(i);
        set(hbg.Nodes(i),'Size',[sz sz]);
    end
end

% Now I can use dolayout:
set(bg,'NodeAutoSize','off','EdgeType','curved',varargin{:});
dolayout(hbg,'pathsOnly',true)