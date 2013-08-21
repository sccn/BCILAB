% [Pdipoles,labels] = label_dipoles(points,pdfgm_regions,goutput)
%
% Description:
%
% Computes the prior probabilities of dipoles to belong to certain
% anatomical regions, given a brain parcelation. 
%
% Inputs: 
%
% points: N X 3 set of N dipole locations in MNI space, if the size of points
% is N X 4, the 4 column represents the localization uncertainty of the
% dipole, this can be used to construct a spherical gaussian around each dipole, by
% default, the sigma of this gaussian is 0 if this value is not provided.  
%
% pdfgm_regions: string list of files containing the maximum probability
% maps of each anatomical region. If this argument is not provided a
% uigetdir box apear allowing you to select the directory where the maps are.
%
% goutput: boolean, 1 (default) => the program produce a graphical output (nice for papers)
%
% Outputs:
%
% Pdipoles: N X Number of anatomical regions, contains the probability of each
% dipole to belong to different anatomical regions.
%
% labels: cell array with the name of the anatomical regions
%
% Note: This software has been tested using the LPBA40/SPM5(GM) maps
% povided by the LONI project and defined in the article: David W. Shattuck, et al.,
% Construction of a 3D probabilistic atlas of human cortical structures, NeuroImage,
% Volume 39, Issue 3, 1 February 2008, Pages 1064-1080, ISSN 1053-8119, DOI: 10.1016/j.neuroimage.2007.09.031.
% (http://www.sciencedirect.com/science/article/B6WNP-4PSC25J-1/2/213a6d22d8a5b67b411e2de09a7d03bf)
%
% Author: Alejandro Ojeda, SCCN, INC, UCSD, 13-Apr-2011
% email: alejandro@sccn.ucsd.edu


function [Pdipoles,labels] = label_dipoles(points,pdfgm_regions,goutput)
error(nargchk(1,3, nargin));
if nargin < 2
    path = [fileparts(which('label_dipoles.m')) filesep 'PDFGM'];
    pdfgm_regions = dir(path);
    pdfgm_regions([1 2]) = [];
    pdfgm_regions = [repmat([path filesep],length(pdfgm_regions),1) char(pdfgm_regions.name)];
    if isempty(pdfgm_regions)
        path = uigetdir(pwd,'Enter the folder where the PDF of the regions are store');
        pdfgm_regions = dir(path);
        pdfgm_regions([1 2]) = [];
        pdfgm_regions = char(pdfgm_regions.name);
        if isempty(pdfgm_regions), error('Invalid PDF folder');end
    end
end
if nargin < 3, goutput = 0;end

%   0 =< bspline_degree =< 7; 
bspline_degree = [1 1 1 0 0 0]; % linear: speed up and avoid noisy interpolation

V = nifti2spm_vol(pdfgm_regions);
Nregions = length(V);

[M,N] = size(points);
if N<3, error('size(points,2) >=3'),end
if N==3, points(:,4) = 0;end
sigma = 7*points(:,4);

Pdipoles = zeros(M,Nregions);
if goutput
    hwait = waitbar(0,'Integrating in dipol''s neighborhood...');
    if logical(sigma), cmd = 'waitbar(it/M,hwait);';else cmd = 'waitbar(jt/Nregions,hwait);';end
else
    cmd = '';
end

if logical(sigma)
    for it=1:M
        
        % generate unit cube
        unidimgrid = -sigma(it):sigma(it);
        [x,y,z] = ndgrid(unidimgrid,unidimgrid,unidimgrid);
        
        % find the points inside the sphere
        R  = sqrt(sum(([x(:) y(:) z(:)]).^2,2));
        I = R<=sigma(it);
        X = [x(:) y(:) z(:)];
        X = X(I,:);
        Np = sum(I);
        dX = ones(Np,1)*points(it,1:3); % ultra fast repmat
        Xi = X + dX;
        Ones = ones(Np,1);
        
        % dipfit prior
        Pdip = normpdf(R(I),0,sigma(it));
        
        % anatomical prior (LONI)
        Preg = zeros(Np,Nregions);
        for jt=1:Nregions
            
            % from MNI to voxel
            Xj = V(jt).mat\[Xi Ones]';
            Xj = Xj(1:3,:)';
            
            % B-spline interpolation in the volume
            Preg(:,jt) = interp_image_volume(V(jt),Xj,bspline_degree);
        end
        
        % eliminating points out of the box
        Preg(isnan(Preg)) = 0; 
        
        % integration in the neighborhood of the point for each region
        Pdipoles(it,:) = sum((Pdip*ones(1,Nregions)).*Preg);
        Pdipoles(it,:) = Pdipoles(it,:)/sum(Pdipoles(it,:));
        eval(cmd);
    end
else
    Ones = ones(M,1);
    Preg = zeros(M,Nregions);
    for jt=1:Nregions
        
        % from MNI to voxel
        Xj = V(jt).mat\[points(:,1:3) Ones]';
        Xj = Xj(1:3,:)';
        
        % B-spline interpolation in the volume
        Preg(:,jt) = interp_image_volume(V(jt),Xj,bspline_degree);
        
        eval(cmd);
    end
    
    % eliminating points out of the box
    Preg(isnan(Preg)) = 0;
    
    % integration in the neighborhood of the point for each region
    Pdipoles = Preg;
    Pdipoles = Pdipoles./(sum(Pdipoles,2)*ones(1,Nregions));
end
Pdipoles(isnan(Pdipoles)) = 0;


labels = cell(Nregions,1);
for jt=1:Nregions
    ind = strfind(V(jt).fname,'.');
    labels{jt} = V(jt).fname(ind(3)+1:ind(end-2)-1);
    ind = strfind(labels{jt},'.');
    labels{jt}(ind) = ' ';
end

% make labels look nicer: capitalize the first letter of each word and replace _  with space
for i=1:length(labels)
    labels{i} = strrep(labels{i}, '_', ' ');
    labels{i} = regexprep(regexprep(labels{i}, '(^.)', '${upper($1)}'), '(?<= \s*)([a-z])','${upper($1)}');
end;

if goutput
    close(hwait);
    figure;
    imagesc(Pdipoles')
    set(gca,'YTick',1:length(labels))
    set(gca,'YTickLabel',labels)
    colorbar
    grid on;
    title('Prior(dipole,anatomical region)')
end
