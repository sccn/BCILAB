function colors = hlp_getROIVertexColorTable(varargin)

% numVerticesInMesh - the number of vertices in the mesh you wish to color
% roiVertices - a cell array where the kth cell contains the vertex indices
%               of the kth ROI.
% bgColor     - a [1 x 3] RGB triplet indicating the color of all vertices
%               that are not found in roiVertices.
% roiColors   - a function name, function handle, or matrix indicating how
%               to color the ROIs defined by roiVertices.
%               If function name or function handle: this function must
%               accept the number of ROIs and return a [num_roi x 3] matrix
%               of RGB triplets.
%               If a matrix, it must be [num_roi x 3] matrix of RGB triplet

arg_define([0 Inf],varargin, ...
    arg_nogui({'numVerticesInMesh','NumVerticesInMesh'},[],[],'The number of vertices in the mesh you wish to color'), ...
    arg_nogui({'roiVertices','RoiVertices'},[],[],'Roi vertices. A cell array where the kth cell contains the vertex indices of the kth ROI.'), ...
    arg({'bgColor','BackgroundColor'},[0.1 0.1 0.1],[],'Background Color. A [1 x 3] RGB triplet indicating the color of all vertices that are not found in roiVertices.'), ...
    arg({'roiColors','RoiColors'},@(x)distinguishable_colors(x,[0.1 0.1 0.1]),[],'Roi Colors. A function name, function handle, or matrix indicating how to color the ROIs defined by roiVertices. If function name or function handle: this function must accept the number of ROIs and return a [num_roi x 3] matrix of RGB triplets. If a matrix, it must be [num_roi x 3] matrix of RGB triplet','type','expression','shape','row'));


nMyROIs = length(roiVertices);

if ischar(roiColors)
    roiColors = feval(roiColors,nMyROIs);
elseif isa(roiColors,'function_handle')
    roiColors = roiColors(nMyROIs);
elseif ismatrix(roiColors) && size(roiColors,1)==1
    roiColors = repmat(roiColors,nMyROIs,1);
end

colors = repmat(bgColor,numVerticesInMesh,1);

for k=1:length(roiVertices)
    colors(roiVertices{k},:) = repmat(roiColors(k,:),length(roiVertices{k}),1);
end