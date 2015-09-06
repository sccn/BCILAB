function dipfit = hlp_makeDipfitStruct(sourceSpace,roiVertices,reducedSpace,CSD,sourceCoords)
% create dipfit structure containing the locations (.posxyz) and moments 
% (.momxyz) of a set of dipoles located at each vertex of the source space.  
% Alternately, one can provide a cell array of integers (roiVertices) where
% each cell contains the indices of vertices within a given region of 
% interest. If this is provided, then dipfit.model(k).posxyz contains the
% centroid location of the kth ROI and momxyz contains its mean moment.
    
% The dipfit structure contains centroids (dipfit.model.posxyz) and
% surface mesh (dipfit.model.surfmesh) for each ROI as well as complete 
% surface mesh (dipfit.surfmesh). We also store the indices of each ROI
% into the complete surface mesh (dipfit.model.meshVertices)
    
% Inputs:
%   sourceSpace:    Surface structure containing fields 
%                     .vertices [num_vertices x 3]
%                     .faces [num_faces x 3]
% Optional:
%   roiVertices:    A cell array of integers (roiVertices) where each cell 
%                   contains the indices of vertices within a given ROI
%   CSD:            [num_vertices x num_samples] current source density
%
% Outputs:
%   dipfit:         dipfit model (struct array) containing fields
%                   .surfmesh:           full surface mesh (faces, vertices)
%                   .model.posxyz:       dipole location (X,Y,Z)
%                   .model.momxyz:       dipole moment
%                   .model.surfmesh:     ROI surface mesh (faces,vertices)
%                   .model.meshVertices: ROI vertex indices into full mesh
%
% Author: Tim Mullen, 2013, SCCN/INC/UCSD

if nargin<1
    error('You must provide a sourceSpace');
end
if nargin<2
    roiVertices = {};
end
if nargin<3
    reducedSpace = [];
end
if nargin<4
    CSD = [];
end
if nargin<5
    sourceCoords = [];
end

% determine number of dipoles
if isempty(roiVertices)
    numDipoles = size(sourceSpace.vertices,1);
else
    numDipoles = length(roiVertices);
end


% compute normal vectors
if ~isempty(CSD)
    normals = geometricTools.getSurfaceNormals(sourceSpace.vertices,sourceSpace.faces,false);
    if size(CSD,1) == 3*size(sourceSpace.vertices,1)  
        CSD = reshape(CSD,[size(CSD,1)/3 3]);
        momxyz  = sqrt(sum(CSD.^2,2));
        s   = sign(dot(normals,CSD,2));
        momxyz  = s.*momxyz; % momxyz are the moment vectors (scaled by power)
    else
        momxyz = CSD;
    end
else
    momxyz = zeros(numDipoles,3);
end

if isempty(sourceCoords)
    posxyz = zeros(numDipoles,3);
else
    posxyz = sourceCoords;
end

% compute ROI averages
if ~isempty(roiVertices)
    
    % average vertex locations (spatial centroid)
    % obtain surfaces for each ROI
    surfmesh = cell(1,numDipoles);
    for k=1:numDipoles
        [v,f]       = geometricTools.getSurfaceROI(sourceSpace.vertices,...
                                                   sourceSpace.faces,   ...
                                                   roiVertices{k});
        surfmesh{k} = struct('vertices',v,'faces',f);
        if isempty(sourceCoords)
            posxyz(k,:) = mean(meshcentroid(v,f));
        end
    end
    
    % average vertex moments
    if ~isempty(CSD)
        momxyz = cellfun(@(x) mean(momxyz(x,:),1)', ...
                     roiVertices, 'UniformOutput',false);
        momxyz = cell2mat(momxyz)';
    end
end
       
% construct the dipole model
dipfit.hdmfile  = '';
dipfit.mrifile  = '';
dipfit.surfmesh = sourceSpace;
dipfit.reducedMesh = reducedSpace;
dipfit.chanfile = '';
dipfit.chansel  = [];
dipfit.coordformat = 'MNI';
dipfit.coord_transform = [0 0 0 0 0 0 1 1 1];

for k=1:numDipoles
    dipfit.model(k).posxyz      = posxyz(k,:);
    dipfit.model(k).momxyz      = momxyz(k,:);
    dipfit.model(k).rv          = 0;
    dipfit.model(k).select      = 1;
    dipfit.model(k).diffmap     = [];
    dipfit.model(k).sourcepot   = [];
    dipfit.model(k).datapot     = [];
    dipfit.model(k).meshVertices= roiVertices{k};
    dipfit.model(k).surfmesh    = surfmesh{k};
end


function centroid=meshcentroid(v,f)
%
% centroid=meshcentroid(v,f)
% 
% compute the centroids of a mesh defined by nodes and elements
% (surface or tetrahedra) in R^n space
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%      v: surface node list, dimension (nn,3)
%      f: surface face element list, dimension (be,3)
%
% output:
%      centroid: centroid positions, one row for each element
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

ec=reshape(v(f(:,1:size(f,2))',:)', [size(v,2) size(f,2) size(f,1)]);
centroid=squeeze(mean(ec,2))';
