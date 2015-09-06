function colors = hlp_getROIColorTable(atlasLabels,roiLabels,bgColor,roiColors)
% atlasLabels - a cell array with (string) labels of all ROIs in the
%               Atlas (in proper order).
% roiLabels   - a cell array with labels of all ROIs you want to color
% bgColor     - a [1 x 3] RGB triplet indicating the color of all Atlas
%               ROIs that are not in the roiLabels list.
% roiColors   - a function name, function handle, or matrix indicating how
%               to color the ROIs that *are* in the roiLabels list.
%               If function name or function handle: this function must
%               accept the number of ROIs and return a [num_roi x 3] matrix
%               of RGB triplets.
%               If a matrix, it must be [num_roi x 3] matrix of RGB triplet

nROIAtlas = size(atlasLabels);
roiInds = find(ismember_bc(atlasLabels,roiLabels));
nMyROIs = length(roiInds);

if ischar(roiColors)
    roiColors = feval(roiColors,nMyROIs);
elseif isa(roiColors,'function_handle')
    roiColors = roiColors(nMyROIs);
elseif ismatrix(roiColors) && size(roiColors,1)==1
    roiColors = repmat(roiColors,nMyROIs,1);
end

colors = repmat(bgColor,nROIAtlas,1);
colors(roiInds,:) = roiColors;


