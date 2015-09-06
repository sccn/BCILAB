function metadata = hlp_nft2mobilab(SubjPath,SubjName,MobilabHeadModelName,chanlabels,atlas,AtlasImgFile,AtlasRoiFile,AtlasNumRoi)
% convert NFT forward model and surfaces to Mobilab object metadata
% This function requires NFT, Mobilab and SPM (for Atlas generation)
%
% Author: Tim Mullen, SCCN/INC 2013

%% set up defaults
if nargin < 3
    error('you must provide at least 3 inputs'); end
if nargin < 4
    chanlabels = []; end
if nargin < 5
    atlas = []; end
if nargin < 6
    AtlasImgFile = '';
elseif nargin < 8
    error('If AtlasImgFile is provided, you must also provide AtlasRoiFile and AtlasNumRoi');
end

RootPath = fullfile(SubjPath,SubjName);

%% load meshes from file
Surfaces = mesh2volstr(RootPath);

%% load LFM
K = load([RootPath,'_LFM.mat']);
K = K.LFM;

%% load sensor locations
sens = load([RootPath, '.sensors'],'-mat');

%% Create Atlas
if ~isempty(AtlasImgFile)
    fprintf('Generating atlas...');
    % map AAL atlas onto this cortical surface
    surfd = struct('faces',Surfaces.bnd(3).tri,'vertices',Surfaces.bnd(3).pnt);
    atlas = geometricTools.labelSurface(surfd,AtlasImgFile, AtlasRoiFile, AtlasNumRoi);
    atlas.color = atlas.colorTable;
    fprintf('done.\n');
end

% create a 'dummy' atlas, if needed
if isempty(atlas)
    atlas = struct('colorTable',ones(size(Surfaces.bnd(3).pnt,1)),'label',{'brain'}); 
end

%% Assemble the head model metadata

if isempty(chanlabels)
    chanlabels = cellstr(num2str((1:size(sens.pnt,1))')); end

if length(chanlabels) ~= size(sens.pnt,1)
    error('number of channel labels must equal number of sensors');
end

% get the fieldnames for a new Mobilab headmodel
metadata = [];
fnames = fieldnames(headModel());
for f=fnames'
    metadata.(f{1}) = [];
end
metadata.atlas          = atlas;
metadata.channelSpace   = sens.pnt;
metadata.label          = chanlabels;
metadata.leadField      = K;
metadata.leadFieldFile  = [MobilabHeadModelName '_LFM.mat'];
metadata.surfacesFilename = [MobilabHeadModelName '_Surfaces.mat'];
metadata.surfaces       = arrayfun(@(bnd) struct('faces',bnd.tri,'vertices',bnd.pnt), Surfaces.bnd);

% % save head model
% save(MobilabHeadModelOutputPath,'metadata');

