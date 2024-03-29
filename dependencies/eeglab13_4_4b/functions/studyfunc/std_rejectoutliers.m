% std_rejectoutliers()  - Commandline function, to reject outlier component(s) from clusters. 
%                            Reassign the outlier component(s) to an outlier cluster specific to each cluster. 
% Usage:    
%                   >> [STUDY] = std_rejectoutliers(STUDY, ALLEEG, clusters, th);   
% Inputs:
%   STUDY         - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG        - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                       ALLEEG for a STUDY set is typically created using load_ALLEEG().  
% Optional inputs:
%   clusters        - [numeric vector| 'all' ] specific cluster numbers (or 'all' clusters), which outliers  
%                        will be rejected from. {default:'all'}.   
%   th                 - [number] a threshold factor to select outliers. How far a component can be from the 
%                       cluster centroid (in the cluster std multiples) befor it will be considered as an outlier. 
%                       Components that their distance from the cluster centroid are more than this factor 
%                       times the cluster std (th *std) will be rejected. {default: 3}.          
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with the components reassignment,
%                    from the cluster to its outlier cluster. 
%
%   Example:
%                         >> clusters = [10 15]; th = 2;   
%                         >> [STUDY] = std_rejectoutliers(STUDY, ALLEEG, clusters, th);  
%                    Reject outlier components (that are more than 2 std from the cluster centroid) from cluster 10 and 15. 
%
%  See also  pop_clustedit         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, July, 2005

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, July 11, 2005, hilit@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function STUDY = std_rejectoutliers(STUDY, ALLEEG, varargin)

cls = 2:length(STUDY.cluster); % all clusters in STUDY
th = 3; % The threshold factor - default: 3 

if length(varargin) > 1
    if isnumeric(varargin{1})
        cls = varargin{1};
        if isempty(cls)
            cls = 2:length(STUDY.cluster);
        end
    else
        if isstr(varargin{1}) & strcmpi(varargin{1}, 'all')
            cls = 2:length(STUDY.cluster);
        else
            error('std_prejectoutliers: clusters input takes either specific clusters (numeric vector) or keyword ''all''.');
        end
    end
end
tmp =[];
for k = 1: length(cls)
    % don't include 'Notclust' clusters
    if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) & ~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)
        tmp = [tmp cls(k)];
    end
end
cls = tmp;
clear tmp

if length(varargin) == 2
    if isnumeric(varargin{2})
        th = varargin{2};
    else
        error('std_prejectoutliers: std input must be a numeric value.');
    end
end

% Perform validity checks
for k = 1:length(cls)
	% Cannot reject outlier components if cluster is a 'Notclust' or 'Outlier' cluster
	if strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) | strncmpi('Outliers',STUDY.cluster(cls(k)).name,8) | ...
        strncmpi('ParentCluster', STUDY.cluster(cls(k)).name,13)
        warndlg2('Cannot reject outlier components from a Notclust or Outliers cluster');
        return;
	end	
	% Cannot reject outlier components if cluster has children clusters
	if ~isempty(STUDY.cluster(cls(k)).child)   
        warndlg2('Cannot reject outlier components if cluster has children clusters.');
        return;
	end
    
    % If the PCA data matrix of the cluster components is empty (case of merged cluster) 
    if isempty(STUDY.cluster(cls(k)).preclust.preclustdata) % No preclustering information
        warndlg2('Cannot reject outlier components if cluster was not a part of pre-clustering.');
        return;
	end
end

% For each of the clusters reject outlier components 
for k = 1:length(cls)
    % The PCA data matrix of the cluster components
    clsPCA = STUDY.cluster(cls(k)).preclust.preclustdata;
    % The cluster centroid
    clsCentr = mean(clsPCA,1);
    % The std of the cluster (based on the distances between all cluster components to the cluster centroid). 
    std_std = std(sum((clsPCA-ones(size(clsPCA,1),1)*clsCentr).^2,2),1);
    outliers = []; 
    for l = 1:length(STUDY.cluster(cls(k)).comps)
        compdist = sum((clsPCA(l,:) - clsCentr).^2); % Component distance from cluster centroid
        if compdist > std_std * th % check if an outlier
            outliers = [ outliers l];
        end
    end
    % Move outlier to the outlier cluster
    if ~isempty(outliers) % reject outliers if exist
        STUDY = std_moveoutlier(STUDY, ALLEEG,cls(k) , outliers); 
    end    
end
