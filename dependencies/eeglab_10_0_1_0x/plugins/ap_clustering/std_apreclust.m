function STUDY = std_apreclust(STUDY,ALLEEG, measuresToUseInClustering, reCalculateAll, correctErpCorrelationWithScalpMapCorrelation)
% std_apreclust() - Calculates pre-clustering measures (pairwise component similarity matrices)
%                   for use in Affinity Product (AF) component clustering.
%                   Results are placed in the field: STUDY.preclust.similarity
% Usage:
%     >> STUDY = std_apreclust(STUDY,ALLEEG, measures, reCalculate);
%
% Inputs:
%
%     STUDY      - STUDY data structure
%     ALLEEG     - ALLEEG vector of loaded EEG data structures
%     measures   - a cell-array of strings telling which measures should be calculated for each STUDY condition:
%
%     'erp' = component ERPs
%     'spec' = component log mean spectra
%     'ersp' = component ERSPs
%     'itc' = component ITCs
%     'dipole' = equivalent dipole locations
%     'map' = component scalp maps ('dipole' may be preferable)
%
% Example:
%
%     >> STUDY = std_apreclust(STUDY,ALLEEG,{'dipole' ,'erp' ,'ersp'});
%     indicates that similarity matrices for component equivalent dipole locations, ERPs, and ERSPs
%     should be pre-computed prior to MP clustering.
%
% Optional Inputs:
%
%     reCalculate - (boolean) If true, the function will remove all previously calculated similarity values
%
%     and re-calculate the similariry matrices.
%
% Output:
%
%     STUDY      - STUDY data structure with STUDY.preclust.similarity field updated
%
%     See also: pop_apreclust(), std_apcluster()
%
% Author: Nima Bigdely-Shamlo and and Brandon Burdge, SCCN/INC/UCSD, 2009


STUDY.preclust.similarity.measuresToUseInClustering = measuresToUseInClustering;

if nargin < 4
    reCalculateAll = true;
end;

if nargin < 5
    correctErpCorrelationWithScalpMapCorrelation = true;
end;

% if isempty(parentClusterID)%isempty(STUDY.etc) || ~isfield(STUDY.etc, 'preclust') || isempty(STUDY.etc.preclust)
%      parentClusterID = 1;
% end;

%if precluster data does not exist, cluster the parent cluster (#1)
if isempty(STUDY.etc) || ~isfield(STUDY.etc, 'preclust') || isempty(STUDY.etc.preclust)
    STUDY.etc.preclust.clustlevel = 1;
end;

% if reCalculateAll
%     STUDY.preclust.similarity = []; % clear all the measure data so they are re-calculated below.
% end;


numberOfConditions = length(STUDY.condition);
numberOfGroups = max(length(STUDY.group), 1); % at least one group is present, even if it is not explicitly defined.
numberOfComponents = length(STUDY.cluster(STUDY.etc.preclust.clustlevel).comps);

% ERP measure
if  ismember('erp', measuresToUseInClustering)
    if (reCalculateAll || ~isfield(STUDY.preclust,'similarity') || ~isfield(STUDY.preclust.similarity,'erpCorr') || size(STUDY.preclust.similarity.erpCorr,1) ~= numberOfComponents)
            
        [STUDY, erpdataConditionByGroup] = std_readerp(STUDY, ALLEEG,  'clusters', 1);        
        
        % perform a 20 Hz lowpass before comparing ERPs, this will reduce the noise and create better
        [b,a] = butter(2,20/(ALLEEG(1).srate / 2),'low');
        
        % combine conditions in each group, then combine groups together
        try
            for group=1:numberOfGroups
                
                % lowpass filter each condition for each gropup
                for i=1:size(erpdataConditionByGroup,1)
                    clustinfo.erpdata{i, group} = filtfilt(b,a, double(erpdataConditionByGroup{i, group}));
                end;
                
                combinedConditions{group} = cat(1,clustinfo.erpdata{:,group});
            end;
        catch exception
            fprintf('Error in concatenating condition. Some conditions for certain subjects might be missing.\n');
            disp(exception.message);
        end;
        
        erpdata = cat(2, combinedConditions{:});
        
        
        % calculate correlations
        fprintf('\nCalculating ERP correlations...');
        erpCorr = 1-squareform(pdist(erpdata','correlation'));                
        
        % change ERP polarity for each pair similarity based on the correlated between their scalpmaps
        % read scalp maps
        if correctErpCorrelationWithScalpMapCorrelation
            fprintf('\nCalculating Scalp maps correlations for correting ERP correlations...');
            grid = zeros(length(STUDY.cluster(1).comps), 3409);
            for i = 1:length(STUDY.cluster(1).comps)
                if mod(i,50) == 0
                    fprintf('%d...',round(100 * i / length(STUDY.cluster(1).comps)));
                end;
                [grid(i,:), y, x ] = std_readtopo(ALLEEG, STUDY.cluster(1).sets(1, i), STUDY.cluster(1).comps(i), 'mode', 'preclust');
            end;
            scalpMapCorrelation =  1-squareform(pdist(grid(:,:),'correlation'));
            
            STUDY.preclust.similarity.erpCorr = erpCorr .* sign(scalpMapCorrelation);
            
            % to calculate euclidean distances, we have to calculate it for each pair separately and
            % consider their scalp-map correlations.
            fprintf('\nCalculating ERP euclidean distances...');
            erpEuclideanDistance = zeros(size(erpCorr));
            for i=1:size(erpEuclideanDistance,1)
                
                if mod(i,50) == 0
                    fprintf('%d...',round(100 * i / length(STUDY.cluster(1).comps)));
                end;
                
                for j=1:(i-1) % only the lower diagonal
                    if sign(scalpMapCorrelation(i,j)) > 0 % if sign is fine (scalp maps correlation is positive)
                        erpEuclideanDistance(i,j) = sqrt(sum((erpdata(:,i) - erpdata(:,j)).^2));
                    else % if sign should be inverted
                        erpEuclideanDistance(i,j) = sqrt(sum((erpdata(:,i) + erpdata(:,j)).^2));
                    end;
                    
                    erpEuclideanDistance(j,i) = erpEuclideanDistance(i,j);
                end; 
            end;
        else 
            STUDY.preclust.similarity.erpCorr = erpCorr;
            fprintf('\nCalculating ERP euclidean distances...');
            erpEuclideanDistance = squareform(pdist(erpdata'));
        end;
        
        STUDY.preclust.similarity.erpEuclideanDistance = erpEuclideanDistance;
    else
        erpCorr =  STUDY.preclust.similarity.erpCorr;
    end;
end;


% ERSP measure
if ismember('ersp', measuresToUseInClustering)
    if (reCalculateAll || ~isfield(STUDY.preclust,'similarity') || ~isfield(STUDY.preclust.similarity,'erspCorr') || size(STUDY.preclust.similarity.erspCorr,1) ~= numberOfComponents)
        fprintf('\n');
            
        [STUDY, erspdataConditionByGroup] = std_readersp(STUDY, ALLEEG,  'clusters', 1);
        
        % combine conditions in each group, then combine groups together
        for group=1:numberOfGroups
            combinedConditions{group} = cat(1, erspdataConditionByGroup{:,group});
        end;
        
        erspdata = cat(3, combinedConditions{:});
        
        % calculate correlations
        fprintf('\nCalculating ERSP correlations and eauclidean distances...');
        
        % swap dimensions and linearize 3D matrix
        erspdataWithCompsFirst = permute(erspdata,[3 2 1]);
        erspdataLinearized = erspdataWithCompsFirst(:,:)';
        
        erspCorr = 1-squareform(pdist(erspdataLinearized','correlation'));
        erspEuclideanDistance = squareform(pdist(erspdataLinearized'));
        
        STUDY.preclust.similarity.erspCorr = erspCorr;
        STUDY.preclust.similarity.erspEuclideanDistance = erspEuclideanDistance;
    else
        erspCorr =  STUDY.preclust.similarity.erspCorr;
    end;
end;



% ITC measure
if ismember('itc', measuresToUseInClustering)
    if (reCalculateAll || ~isfield(STUDY.preclust,'similarity') || ~isfield(STUDY.preclust.similarity,'itcCorr') || size(STUDY.preclust.similarity.itcCorr,1) ~= numberOfComponents)
        fprintf('\n');
        
        %         [STUDY clustinfo] = std_readdata(STUDY,ALLEEG,'clusters',STUDY.etc.preclust.clustlevel,'infotype','itc');
        %         clustinfo = STUDY.cluster(STUDY.etc.preclust.clustlevel);
        
        try
            [STUDY, itcdataConditionByGroup] = std_readitc(STUDY, ALLEEG,  'clusters', 1);
        catch
            [STUDY itcdataConditionByGroup] = std_readersp(STUDY, ALLEEG, 'clusters',1, 'infotype','itc');
        end;
        
        % combine conditions in each group, then combine groups together
        for group=1:numberOfGroups
            combinedConditions{group} = cat(1, itcdataConditionByGroup{:,group});
        end;
        itcdata = cat(3, combinedConditions{:});
        
        % calculate correlations
        fprintf('\nCalculating ITC correlations and eauclidean distances...');
        
        % swap dimensions and linearize 3D matrix
        itcdataWithCompsFirst = permute(itcdata,[3 2 1]);
        itcdataLinearized = itcdataWithCompsFirst(:,:)';
        
        itcCorr = 1-squareform(pdist(itcdataLinearized','correlation'));
        itcEuclideanDistance = squareform(pdist(itcdataLinearized'));
        
        STUDY.preclust.similarity.itcCorr = itcCorr;
        STUDY.preclust.similarity.itcEuclideanDistance = itcEuclideanDistance;
    else
        itcCorr =  STUDY.preclust.similarity.itcCorr;
    end;
end;


% DIPOLE measure
if ismember('dipole', measuresToUseInClustering)
   if (reCalculateAll || ~isfield(STUDY.preclust,'similarity') || ~isfield(STUDY.preclust.similarity,'compDistance') || size(STUDY.preclust.similarity.compDistance,1) ~= numberOfComponents)
      fprintf('\n');
      [STUDY clustinfo] = std_readdata(STUDY,ALLEEG,'clusters',STUDY.etc.preclust.clustlevel,'infotype','dipole');
      clustinfo = STUDY.cluster(STUDY.etc.preclust.clustlevel);

      if isfield(clustinfo, 'alldipoles')
         dipoles = clustinfo.alldipoles;
      else
         % combine conditions in each group, then combine groups together
         for group=1:numberOfGroups
            %            combinedConditions{group} = cat(1,clustinfo.dipoles{:,group});
            combinedConditions{group} = cat(1,clustinfo.dipoles{1,group}); % only dipoles from the first condition
         end;
         dipoles = cat(2, combinedConditions{:});
      end;


        % calculate distance between dipoles
        compPos = [];compDistance = [];
        
        for i=1:size(dipoles,2)
            for j=1:size(dipoles,2)
                compDistance(i,j) = norm(dipoles(i).posxyz(1,:)- dipoles(j).posxyz(1,:),'fro');
                %          compDistance(i,j) = distanceBetweenDipoles(dipoles(i).posxyz, dipoles(j).posxyz);
            end;
            compPos(i,:) = dipoles(i).posxyz(1,:);
        end;
        fprintf('\nDipoles done.');
        STUDY.preclust.similarity.compPos = compPos;
        STUDY.preclust.similarity.compDistance = compDistance;
    else
        compPos = STUDY.preclust.similarity.compPos;
        compDistance = STUDY.preclust.similarity.compDistance;
    end;
end;


% Spectra measure
if ismember('spec', measuresToUseInClustering)
    if (reCalculateAll || ~isfield(STUDY.preclust,'similarity') || ~isfield(STUDY.preclust.similarity,'specCorr') || size(STUDY.preclust.similarity.specCorr,1) ~= numberOfComponents)
        fprintf('\nSpectra: percent done = ');
        
        try
           [STUDY specdata ] = std_readspec(STUDY, ALLEEG, 'clusters', 1);
            specdata = cat(1, specdata{:});
        catch
           [STUDY clustinfo] = std_readdata(STUDY,ALLEEG,'clusters',STUDY.etc.preclust.clustlevel,'infotype','spec');
           clustinfo = STUDY.cluster(STUDY.etc.preclust.clustlevel);
           %   specdata = cat(2,clustinfo.specdata{1}, clustinfo.specdata{2});

           % combine conditions in each group, then combine groups together
           for group=1:numberOfGroups
              combinedConditions{group} = cat(1,clustinfo.specdata{:,group});
           end;
           specdata = cat(2, combinedConditions{:});
        end;
        
        specdataEachMeanRemoved = specdata - repmat(mean(specdata, 1), size(specdata,1), 1);
        meanSpec= mean(specdataEachMeanRemoved, 2);
        normalizedSpec = specdataEachMeanRemoved - repmat(meanSpec, 1, size(specdata,2));
        
        
        % calculate correlations
        fprintf('\nCalculating Spectrum correlations...');
        specCorr = 1-squareform(pdist(specdata','correlation'));
        normalizedSpecCorr = 1-squareform(pdist(normalizedSpec','correlation'));
        
        STUDY.preclust.similarity.specCorr = specCorr;
        STUDY.preclust.similarity.normalizedSpecCorr = normalizedSpecCorr;
        STUDY.preclust.similarity.normalizedSpec = normalizedSpec;
        
    else
        specCorr =  STUDY.preclust.similarity.specCorr;
        normalizedSpecCorr = STUDY.preclust.similarity.normalizedSpecCorr;
    end;
end;




% Scalp measure
if ismember('map', measuresToUseInClustering)
    if (reCalculateAll || ~isfield(STUDY.preclust,'similarity') || ~isfield(STUDY.preclust.similarity,'mapCorr') || size(STUDY.preclust.similarity.mapCorr,1) ~= numberOfComponents)
        
        [STUDY clustinfo] = std_readdata(STUDY,ALLEEG,'clusters',STUDY.etc.preclust.clustlevel,'infotype','map');
        clustinfo = STUDY.cluster(STUDY.etc.preclust.clustlevel);
        fprintf('\nScalp maps: percent done = ');
        
        mapdata = zeros(numel(clustinfo.topoall{1}),numberOfComponents );
        for i =1:numberOfComponents
            mapdata(:,i) = clustinfo.topoall{i}(:);
        end;
        
        % removing NAN indices
        s = sum(mapdata,2);
        mapdata(find(isnan(s)),:) = [];
        
        mapCorr = [];
        for i=1:size(mapdata,2)
            if mod(i, 50)  == 0
                fprintf('%d...',round(100*i/size(mapdata,2)));
            end;
            
            for j=1:size(mapdata,2)
                c = corrcoef(mapdata(:,i), mapdata(:,j));
                mapCorr(i,j) = c(1,2);
            end;
        end;
        STUDY.preclust.similarity.mapCorr = mapCorr;
    else
        mapCorr =  STUDY.preclust.similarity.mapCorr;
    end;
end;

fprintf('\n');




cluster_ind = STUDY.etc.preclust.clustlevel;

% of interest for each set of condition
for k = 1:size(STUDY.setind,2)
    % Find the first entry in STUDY.setind(:,k) that is non-NaN. We only need one since
    % components are the same across conditions.
    for ri = 1:size(STUDY.setind,1)
        if ~isnan(STUDY.setind(ri,k)), break; end
    end
    sind = find(STUDY.cluster(cluster_ind).sets(ri,:) == STUDY.setind(ri,k));
    succompind{k} = STUDY.cluster(cluster_ind).comps(sind);
end;
for ind = 1:size(STUDY.setind,2)
    succompind{ind} = succompind{ind}(find(succompind{ind})); % remove zeros
    % (though there should not be any? -Arno)
    succompind{ind} = sort(succompind{ind}); % sort the components
end;

STUDY.etc.preclust.preclustcomps = succompind;
STUDY.etc.preclust.preclustparams = [];
STUDY.etc.preclust.preclustdata = [];

end

function result = distanceBetweenDipoles(posA,posB)
if size(posA,1) == 1 && size(posB,1)==3
    result = norm(posA-posB,'fro');
else
    if size(posA,1) == 1
        posA = repmat(posA,2,1);
    end;
    
    if size(posB,1) == 1
        posB = repmat(posB,2,1);
    end;
    
    d = squareform(pdist([posA;posB]));
    result = min([d(1,3) d(1,4) d(2,3) d(2,4)]); % minimum of pairwise dipoles no belonging to the same IC
end
end
