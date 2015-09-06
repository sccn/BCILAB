% Defines the class headModel for solving forward/inverse problem of the EEG. 
% This class is part of MoBILAB software. 
% For more details visit:  https://code.google.com/p/mobilab/
% 
% Author: Alejandro Ojeda, SCCN/INC/UCSD, Jan-2012

classdef headModel < handle
    properties(GetAccess=public, SetAccess=public,SetObservable)
        channelSpace = [];     % xyz coordinates of the sensors.
        
        fiducials = [];        % xyz of the fiducial landmarks: nassion, lpa, rpa, vertex, and inion.
        
        surfaces = [];         % Pointer to the file where the surfaces representing different layers 
                               % of tissue are stored. The surfaces must be in an array of MATLAB patches
                               % in the following order: 1) scalp, 2) skull, 3) brain (gray matter or 
                               % average between gray and white matter)
                               
        atlas                  % Atlas that labels each vertex in the most internal surface (gray matter).
        
        leadFieldFile = [];    % Pointer to the file where the lead field matrix was stored.
    end
    properties(GetAccess = private, SetAccess = private)
        label;
        F = [];
    end
    properties(Dependent)
        channelLabel = []
    end
    methods
        function obj = headModel(varargin)
            if length(varargin)==1, varargin = varargin{1};end
            
            if ~iscell(varargin) && ischar(varargin) && exist(varargin,'file')
                [obj.channelSpace,obj.label,obj.fiducials] = readMontage(varargin);
                return
            end
            
            ind = find(ismember(varargin(1:2:length(varargin)-1),'channelSpace'));
            if ~isempty(ind), obj.channelSpace = varargin{ind*2};end
            
            ind = find(ismember(varargin(1:2:length(varargin)-1),'surfaces'));
            if ~isempty(ind), obj.surfaces = varargin{ind*2};end
            if ~isempty(obj.surfaces)
                [~,~,e] = fileparts(obj.surfaces);
                if isempty(e), obj.surfaces = [obj.surfaces,'.mat'];end
            end            
            ind = find(ismember(varargin(1:2:length(varargin)-1),'atlas'));
            if ~isempty(ind)
                tmpAtlas = varargin{ind*2};
                if isfield(tmpAtlas,'color'),
                    tmpAtlas.colorTable = tmpAtlas.color;
                    tmpAtlas = rmfield(tmpAtlas,'color');
                end
                obj.atlas = tmpAtlas;
            end
            ind = find(ismember(varargin(1:2:length(varargin)-1),'fiducials'));
            if ~isempty(ind), obj.fiducials = varargin{ind*2};end
            
            ind = find(ismember(varargin(1:2:length(varargin)-1),'leadFieldFile'));
            if ~isempty(ind), obj.leadFieldFile = varargin{ind*2};end
            ind = find(ismember(varargin(1:2:length(varargin)-1),'label'));
            if ~isempty(ind),
                obj.label = varargin{ind*2};
            else
                N = size(obj.channelSpace,1);
                labels = num2str((1:N)');
                obj.label =num2cell(labels',[N,1])';
                for it=1:N, obj.label{it} = deblank(obj.label{it});end
            end
        end
        function labels = getChannelLabels(obj)
            labels = obj.label;
            warning('This method will be deprecated in the future, instead you can access directly the property channelLabel.')
        end
        function channelLabel = get.channelLabel(obj)
            channelLabel = obj.label;
        end
        %%
        function [roiname,roinumber] = labelDipole(obj,dipole)
            if isempty(obj.F)
                load(obj.surfaces)
                obj.F = scatteredInterpolant(surfData(end).vertices(:,1),...
                    surfData(end).vertices(:,2),surfData(end).vertices(:,3),...
                    obj.atlas.colorTable,'nearest');
            end
            roinumber = obj.F(dipole(:,1),dipole(:,2),dipole(:,3));
            roiname = obj.atlas.label(roinumber);
        end
        %%
        function chanlocs = makeChanlocs(obj)
            % make EEGLAB chanlocs structure from channel locations and
            % labels
            if isempty(which('convertlocs'))
                error('EEGLAB function convertlocs.m is missing.');
            end
            for k=1:length(obj.label)
                chanlocs(k) = struct('labels',obj.label{k}, ...
                                     'ref','', ...
                                     'theta',[], ...
                                     'radius',[], ...
                                     'X',obj.channelSpace(k,1), ...
                                     'Y',obj.channelSpace(k,2), ...
                                     'Z',obj.channelSpace(k,3), ...
                                     'sph_theta', [], ...
                                     'sph_phi',[], ...
                                     'sph_radius',[], ...
                                     'type', 'EEG', ...
                                     'urchan', []);
            end
            chanlocs = convertlocs( chanlocs, 'cart2all');
        end
        %%
        function h = plotHeadModel(obj,~) % do not remove the circumflex, I'm passging a second arguments when this method is called from MoBILAB's gui
            % Plots the different layers of tissue, the sensor positions, and their labels.
            % It colors different regions of the cortical surface according to a defined
            % anatomical atlas. Several interactive options for customizing the figure are
            % available.
            
            if isempty(obj.channelSpace) || isempty(obj.label) || isempty(obj.surfaces);
                error('Head model is incomplete or missing.');
            end
            h = headModelViewerHandle(obj,obj.label);
        end
       %%
        function h = plotMontage(obj,showNewfig)
            % Plots a figure with the xyz distribution of sensors, fiducial landmarks, and
            % coordinate axes.
            
            if isempty(obj.channelSpace) || isempty(obj.label);error('MoBILAB:noChannelSpace','Channel space is empty.');end
            if nargin < 2, showNewfig = true;end
            
            if isa(obj,'eeg')
                color = [0.93 0.96 1];
            else
                color = [0.76 0.77 1];
            end
            if showNewfig, figure('Color',color);end
            h = scatter3(obj.channelSpace(:,1),obj.channelSpace(:,2),obj.channelSpace(:,3),'filled',...
                'MarkerEdgeColor','k','MarkerFaceColor','y','parent',gca);
            
            hold on;
            N = length(obj.label);
            k = 1.1;
            for it=1:N, text('Position',k*obj.channelSpace(it,:),'String',obj.label{it});end
            mx = max(obj.channelSpace);
            k = 1.2;
            line([0 k*mx(1)],[0 0],[0 0],'LineStyle','-.','Color','b','LineWidth',2)
            line([0 0],[0 k*mx(2)],[0 0],'LineStyle','-.','Color','g','LineWidth',2)
            line([0 0],[0 0],[0 k*mx(3)],'LineStyle','-.','Color','r','LineWidth',2)
            text('Position',[k*mx(1) 0 0],'String','X','FontSize',12,'FontWeight','bold','Color','b')
            text('Position',[0 k*mx(2) 0],'String','Y','FontSize',12,'FontWeight','bold','Color','g')
            text('Position',[0 0 k*mx(3)],'String','Z','FontSize',12,'FontWeight','bold','Color','r')
            
            try %#ok
                scatter3(obj.fiducials.nasion(1),obj.fiducials.nasion(2),obj.fiducials.nasion(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.nasion,'String','Nas','FontSize',12,'FontWeight','bold','Color','k');
                
                scatter3(obj.fiducials.lpa(1),obj.fiducials.lpa(2),obj.fiducials.lpa(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.lpa,'String','LPA','FontSize',12,'FontWeight','bold','Color','k');
                
                scatter3(obj.fiducials.rpa(1),obj.fiducials.rpa(2),obj.fiducials.rpa(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.rpa,'String','RPA','FontSize',12,'FontWeight','bold','Color','k');
                
                scatter3(obj.fiducials.vertex(1),obj.fiducials.vertex(2),obj.fiducials.vertex(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.vertex,'String','Ver','FontSize',12,'FontWeight','bold','Color','k');
                
                scatter3(obj.fiducials.inion(1),obj.fiducials.inion(2),obj.fiducials.inion(3),'filled','MarkerEdgeColor','k','MarkerFaceColor','K');
                text('Position',1.1*obj.fiducials.inion,'String','Ini','FontSize',12,'FontWeight','bold','Color','k');
            end
            
            % box on;
            hold off;
            axis equal
            axis vis3d
            grid on;
        end
       %%
        function individualHeadModelFile = warpTemplate2channelSpace(obj,headModelFile,individualHeadModelFile)
            % Warps a template head model to the space defined by the sensor positions (channelSpace). It uses Dirk-Jan Kroon's
            % nonrigid_version23 toolbox.
            %
            % For more details see: http://www.mathworks.com/matlabcentral/fileexchange/20057-b-spline-grid-image-and-point-based-registration
            % 
            % Input arguments:
            %       headModelFile:           pointer to the template head model file. To see an example of
            %                                templates see the folder mobilab/data/headModelXX.mat
            %       individualHeadModelFile: pointer to the warped head model (output file)
            % 
            % Output arguments:
            %       individualHeadModelFile: pointer to the warped head model (same as the second input argument)
            %
            % References: 
            %    D. Rueckert et al. "Nonrigid Registration Using Free-Form Deformations: Application to Breast MR Images".
            %    Seungyong Lee, George Wolberg, and Sung Yong Shing, "Scattered Data interpolation with Multilevel B-splines"

            if nargin < 2, error('Reference head model is missing.');end
            if nargin < 3, individualHeadModelFile = [tempname '.mat'];end
            if isempty(obj.channelSpace) || isempty(obj.label), error('Channel space or labels are missing.');end
            if ~exist(headModelFile,'file'), error('The file you''ve entered does not exist.');end
            
            template = load(headModelFile);
            gTools = geometricTools;
            th = norminv(0.90);
            % mapping source to target spaces: S->T
            % target space: individual geometry
            
            try
                T = [obj.fiducials.nasion;...
                    obj.fiducials.lpa;...
                    obj.fiducials.rpa];
                
                % source space: template
                S = [template.fiducials.nasion;...
                    template.fiducials.lpa;...
                    template.fiducials.rpa;...
                    template.fiducials.vertex];
                
                % estimates vertex if is missing
                if isfield(obj.fiducials,'vertex')
                    if numel(obj.fiducials.vertex) == 3
                        T = [T;obj.fiducials.vertex];
                    else
                        point = 0.5*(obj.fiducials.lpa + obj.fiducials.rpa);
                        point = ones(50,1)*point;
                        point(:,3) = linspace(point(3),1.5*max(obj.channelSpace(:,3)),50)';
                        [~,d] = gTools.nearestNeighbor(obj.channelSpace,point);
                        [~,loc] = min(d);
                        point = point(loc,:);
                        T = [T;point];
                    end
                else
                    point = 0.5*(obj.fiducials.lpa + obj.fiducials.rpa);
                    point = ones(50,1)*point;
                    point(:,3) = linspace(point(3),1.5*max(obj.channelSpace(:,3)),50)';
                    [~,d] = gTools.nearestNeighbor(obj.channelSpace,point);
                    [~,loc] = min(d);
                    point = point(loc,:);
                    T = [T;point];
                end
                
                if isfield(obj.fiducials,'inion')
                    if numel(obj.fiducials.vertex) == 3
                        T = [T;obj.fiducials.inion];
                        S = [S;template.fiducials.inion];
                    end
                end
            catch
                disp('Fiducials are missing in the individual head model, selecting the common set of points based on the channel labels.')
                [~,loc1,loc2] = intersect(obj.getChannelLabels,template.label,'stable');
                T = obj.channelSpace(loc1,:);
                S = template.channelSpace(loc2,:);
            end
            try obj.initStatusbar(1,8,'Co-registering...');end %#ok
            
            % affine co-registration
            [Aff,~,scale] = gTools.affineMapping(S,T);
            if isa(obj,'eeg'), obj.statusbar(1);end
            
            % b-spline co-registration (only fiducial landmarks)
            options.Verbose = true;
            options.MaxRef = 2;
            surfData = template.surfData;
            Ns = length(surfData);
            for it=1:Ns
                surfData(it).vertices = gTools.applyAffineMapping(template.surfData(it).vertices,Aff);
            end
            Saff = gTools.applyAffineMapping(S,Aff);
            [Def,spacing,offset] = gTools.bSplineMapping(Saff,T,surfData(1).vertices,options);
            try obj.statusbar(2);end %#ok
            
            % b-spline co-registration (second pass)
            for it=1:Ns
                surfData(it).vertices = gTools.applyBSplineMapping(Def,spacing,offset,surfData(it).vertices);
            end
            T = obj.channelSpace;
            T(T(:,3) <= min(surfData(1).vertices(:,3)),:) = [];
            [S,d] = gTools.nearestNeighbor(surfData(1).vertices,T);
            z = zscore(d);
            S(abs(z)>th,:) = [];
            T(abs(z)>th,:) = [];
            [Def,spacing,offset] = gTools.bSplineMapping(S,T,surfData(1).vertices,options);
            try obj.statusbar(3);end %#ok
            
            % b-spline co-registration (third pass)
            for it=1:Ns
                surfData(it).vertices = gTools.applyBSplineMapping(Def,spacing,offset,surfData(it).vertices);
            end
            T = obj.channelSpace;
            T(T(:,3) <= min(surfData(1).vertices(:,3)),:) = [];
            [S,d] = gTools.nearestNeighbor(surfData(1).vertices,T);
            z = zscore(d);
            S(abs(z)>th,:) = [];
            T(abs(z)>th,:) = [];
            Tm = 0.5*(T+S);
            [Def,spacing,offset] = gTools.bSplineMapping(S,Tm,surfData(1).vertices,options);
            try obj.statusbar(4);end %#ok
            
            % apply the final transformation
            for it=1:Ns
                surfData(it).vertices = gTools.applyBSplineMapping(Def,spacing,offset,surfData(it).vertices);
                surfData(it).vertices = gTools.smoothSurface(surfData(it).vertices,surfData(it).faces);
            end
            
            % fixing topological defects
            try obj.container.container.statusBar.setText('Fixing topological defects...');end %#ok
            dmax = ones(Ns-1,1)*5;
            dmax(1) = 8;
            dmax = dmax*scale;
            ind = fliplr(1:Ns);
            for it=1:Ns-1
                surfData(ind(it+1)).vertices = gTools.repareIntersectedSurface(surfData(ind(it)),surfData(ind(it+1)),dmax(it));
                try obj.statusbar(it+5);end %#ok
            end
            
            ind =  obj.channelSpace(:,3) > min(surfData(1).vertices(:,3));
            T = gTools.nearestNeighbor(surfData(1).vertices,obj.channelSpace);
            channelSpace = obj.channelSpace; %#ok
            channelSpace(ind,:) = T(ind,:);  %#ok
            [~,loc] = unique(channelSpace,'rows');%#ok
            indInterp = setdiff(1:size(obj.channelSpace,1),loc);
            if ~isempty(indInterp)
                x = setdiff(channelSpace,channelSpace(indInterp,:),'rows');%#ok
                xi = gTools.nearestNeighbor(x,channelSpace(indInterp,:));%#ok
                channelSpace(indInterp,:) = 0.5*(xi + channelSpace(indInterp,:));%#ok
            end
            obj.channelSpace = channelSpace; %#ok
            
            if isfield(template,'atlas'), 
                if isfield(template.atlas,'color')
                    colorTable = template.atlas.color;
                    template.atlas = rmfield(template.atlas,'color');
                    template.atlas.colorTable = colorTable;
                end
                obj.atlas = template.atlas;
            end
            if exist(obj.surfaces,'file'), delete(obj.surfaces);end
            obj.surfaces = individualHeadModelFile;
            save(obj.surfaces,'surfData');
            try obj.statusbar(8);end %#ok
            disp('Done!')
        end
       %%
        function hFigureObj = plotOnModel(obj,J,V,figureTitle)
            % Plots cortical/topographical maps onto the cortical/scalp surface.
            % 
            % Input parameters:
            %       J:           cortical map size number of vertices of the cortical surface by number of time points
            %       V:           topographic map size number of vertices of the scalp surface by number of time points; 
            %                    if V is empty, a single color is used simulating the color of the skin 
            %       figureTitle: title of the figure (optional)
            %                    
            % Output argument:   
            %       hFigure:     figure handle 
            
            if nargin < 2, error('Not enough input arguments');end
            if nargin < 3, V = [];end
            if nargin < 4, figureTitle = '';end
            if isa(obj,'pcdStream'), channelLabels = obj.parent.label;else channelLabels = obj.label;end
            hFigureObj = currentSourceViewer(obj,J,V,figureTitle,channelLabels);
        end
       %%
        function Aff = warpChannelSpace2Template(obj,headModelFile,individualHeadModelFile,regType)
            % Estimates a mapping from channel space to a template's head. It uses Dirk-Jan Kroon's
            % nonrigid_version23 toolbox.
            %
            % For more details see: http://www.mathworks.com/matlabcentral/fileexchange/20057-b-spline-grid-image-and-point-based-registration
            %
            % Input arguments:
            %       headModelFile:           pointer to the template head model file. To see an example
            %                                of templates see the folder mobilab/data/headModelXX.mat
            %       individualHeadModelFile: pointer to the warped head model (output file)
            %       regType:                 co-registration type, could be 'affine' or 'bspline'. In case
            %                                of 'affine' only the affine mapping is estimated (rotation,
            %                                traslation, and scaling). 'bspline' starts from the affine 
            %                                mapping and goes on to estimate a non-linear defformation
            %                                field that captures better the shape of the head.
            %
            % Output arguments:
            %       Aff: affine matrix
            %
            % References: 
            %    D. Rueckert et al. "Nonrigid Registration Using Free-Form Deformations: Application to Breast MR Images".
            %    Seungyong Lee, George Wolberg, and Sung Yong Shing, "Scattered Data interpolation with Multilevel B-splines"

            if nargin < 2, error('Reference head model is missing.');end
            if nargin < 3, individualHeadModelFile = ['surfaces_' num2str(round(1e5*rand)) '.mat'];end
            if nargin < 4, regType = 'bspline';end
            if isempty(obj.channelSpace) || isempty(obj.label) || isempty(obj.fiducials), error('Channel space or fiducials are missing.');end
            if ~exist(headModelFile,'file'), error('The file you''ve entered does not exist.');end
                       
            template = load(headModelFile);
            surfData = template.surfData;
            gTools = geometricTools;
            th = norminv(0.90);
            % mapping source to target spaces: S->T
            % target space: template
            T = [template.fiducials.nasion;...
                template.fiducials.lpa;...
                template.fiducials.rpa;...
                template.fiducials.vertex];
            
            % source space: individual geometry
            S = [obj.fiducials.nasion;...
                obj.fiducials.lpa;...
                obj.fiducials.rpa];
            
            % estimates vertex if is missing
            if isfield(obj.fiducials,'vertex')
                if numel(obj.fiducials.vertex) == 3
                    S = [S;obj.fiducials.vertex];
                else
                    point = 0.5*(obj.fiducials.lpa + obj.fiducials.rpa);
                    point = ones(50,1)*point;
                    point(:,3) = linspace(point(3),1.5*max(obj.channelSpace(:,3)),50)';
                    [~,d] = gTools.nearestNeighbor(obj.channelSpace,point);
                    [~,loc] = min(d);
                    point = point(loc,:);
                    S = [S;point];
                end
            else
                point = 0.5*(obj.fiducials.lpa + obj.fiducials.rpa);
                point = ones(50,1)*point;
                point(:,3) = linspace(point(3),1.5*max(obj.channelSpace(:,3)),50)';
                [~,d] = gTools.nearestNeighbor(obj.channelSpace,point);
                [~,loc] = min(d);
                point = point(loc,:);
                S = [S;point];
            end
            
            if isfield(obj.fiducials,'inion')
                if numel(obj.fiducials.vertex) == 3
                    S = [S;obj.fiducials.inion];
                    T = [T;template.fiducials.inion];
                end
            end
            if isa(obj,'eeg')
                obj.initStatusbar(1,8,'Co-registering...');
            else
                disp('Co-registering...');
            end
            
            % affine co-registration
            Aff = gTools.affineMapping(S,T);
            if isa(obj,'eeg'), obj.statusbar(1);end
            
            obj.channelSpace = gTools.applyAffineMapping(obj.channelSpace,Aff);
            obj.fiducials.lpa = gTools.applyAffineMapping(obj.fiducials.lpa,Aff);
            obj.fiducials.rpa = gTools.applyAffineMapping(obj.fiducials.rpa,Aff);
            obj.fiducials.nasion = gTools.applyAffineMapping(obj.fiducials.nasion,Aff);
            
            if ~strcmp(regType,'affine')
                % b-spline co-registration (only fiducial landmarks)
                options.Verbose = true;
                options.MaxRef = 2;
                Saff = gTools.applyAffineMapping(S,Aff);
                [Def,spacing,offset] = gTools.bSplineMapping(Saff,T,obj.channelSpace,options);
                if isa(obj,'eeg'), obj.statusbar(2);end
                
                % b-spline co-registration (second pass)
                obj.channelSpace = gTools.applyBSplineMapping(Def,spacing,offset,obj.channelSpace);
                obj.fiducials.lpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.lpa);
                obj.fiducials.rpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.rpa);
                obj.fiducials.nasion = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.nasion);
                
                T = template.surfData(1).vertices;
                S = obj.channelSpace;
                S(S(:,3) <= min(T(:,3)),:) = [];
                [S,d] = gTools.nearestNeighbor(S,T);
                z = zscore(d);
                S(abs(z)>th,:) = [];
                T(abs(z)>th,:) = [];
                [Def,spacing,offset] = gTools.bSplineMapping(S,T,obj.channelSpace,options);
                if isa(obj,'eeg'), obj.statusbar(3);end
                
                % b-spline co-registration (third pass)
                obj.channelSpace = gTools.applyBSplineMapping(Def,spacing,offset,obj.channelSpace);
                obj.fiducials.lpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.lpa);
                obj.fiducials.rpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.rpa);
                obj.fiducials.nasion = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.nasion);
                
                T = template.surfData(1).vertices;
                S = obj.channelSpace;
                S(S(:,3) <= min(T(:,3)),:) = [];
                [S,d] = gTools.nearestNeighbor(S,T);
                z = zscore(d);
                S(abs(z)>th,:) = [];
                T(abs(z)>th,:) = [];
                Tm = 0.5*(T+S);
                [Def,spacing,offset] = gTools.bSplineMapping(S,Tm,obj.channelSpace,options);
                if isa(obj,'eeg'), obj.statusbar(4);end
                
                % apply the final transformation
                obj.channelSpace = gTools.applyBSplineMapping(Def,spacing,offset,obj.channelSpace);
                obj.fiducials.lpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.lpa);
                obj.fiducials.rpa = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.rpa);
                obj.fiducials.nasion = gTools.applyBSplineMapping(Def,spacing,offset,obj.fiducials.nasion);
            end
            
            % fixing topological defects
            if isa(obj,'eeg')
                obj.statusbar.setText('Fixing topological defects...');
            else
                disp('Fixing topological defects...');
            end
            Ns = length(surfData);
            dmax = ones(Ns-1,1)*5;
            dmax(1) = 8;
            ind = fliplr(1:Ns);
            for it=1:Ns-1
                surfData(ind(it+1)).vertices = gTools.repareIntersectedSurface(surfData(ind(it)),surfData(ind(it+1)),dmax(it));
                if isa(obj,'eeg'), obj.statusbar(it+5);end
            end
            
            ind =  obj.channelSpace(:,3) > min(surfData(1).vertices(:,3));
            T = gTools.nearestNeighbor(surfData(1).vertices,obj.channelSpace);
            channelSpace = obj.channelSpace; %#ok
            channelSpace(ind,:) = T(ind,:);  %#ok
            [~,loc] = unique(channelSpace,'rows');%#ok
            indInterp = setdiff(1:size(obj.channelSpace,1),loc);
            if ~isempty(indInterp)
                x = setdiff(channelSpace,channelSpace(indInterp,:),'rows');%#ok
                xi = gTools.nearestNeighbor(x,channelSpace(indInterp,:));%#ok
                channelSpace(indInterp,:) = 0.5*(xi + channelSpace(indInterp,:));%#ok
            end
            obj.channelSpace = channelSpace; %#ok
            
            if isfield(template,'atlas'), obj.atlas = template.atlas;end
            if exist(obj.surfaces,'file'), delete(obj.surfaces);end
            obj.surfaces = individualHeadModelFile;
            save(obj.surfaces,'surfData');
            if isa(obj,'eeg'), obj.statusbar(8);end
        end
       %%
        function computeLeadFieldBEM(obj, conductivity,orientation)
            % Computes the lead field matrix interfacing OpenMEEG toolbox [1].
            %
            % Input arguments:
            %       conductivity: conductivity of each layer of tissue, scalp - skull - brain,
            %                     default: 0.33-0.022-0.33 S/m. See [2, 3, 4] for details.
            %        orientation: if true, computes the orientation free lead field, otherwise
            %                     it constrain the dipoles to be normal to the cortical surface
            %
            % The computed lead field is stored inside the object in obj.leadFieldFile.
            %
            % References:
            %   [1] Gramfort, A., Papadopoulo, T., Olivi, E., & Clerc, M. (2010).
            %         OpenMEEG: opensource software for quasistatic bioelectromagnetics.
            %         Biomedical engineering online, 9, 45. doi:10.1186/1475-925X-9-45
            %   [2] Vald??s-Hern??ndez, P.A., Von Ellenrieder, N., Ojeda-Gonzalez, A., Kochen, S.,
            %         Alem??n-G??mez, Y., Muravchik, C., & A Vald??s-Sosa, P. (2009). Approximate
            %         average head models for EEG source imaging. Journal of Neuroscience Methods,
            %         185(1), 125???132.
            %   [3] Wendel, K., Malmivuo, J., 2006. Correlation between live and post mortem skull
            %         conductivity measurements. Conf Proc IEEE Eng Med Biol Soc 1, 4285-4288.
            %   [4] Oostendorp, T.F., Delbeke, J., Stegeman, D.F., 2000. The conductivity of the 
            %         human skull: Results of in vivo and in vitro measurements. Ieee Transactions
            %         on Biomedical Engineering 47, 1487-1492.
                        
            dispCommand = false;
            if nargin < 2, conductivity = [0.33 0.022 0.33];end
            if nargin < 3, orientation = true;end
            if isempty(obj.channelSpace) || isempty(obj.label) || isempty(obj.surfaces);
                error('Head model is incomplete or missing.');
            end
            if any(conductivity == -1)
                prefObj = [...
                    PropertyGridField('conductivity',[0.33 0.022 0.33],'DisplayName','Conductivity','Description',sprintf('Conductivity values are taken from Valdes-Hernandez et al., 2006, check \nalso Oostendrop TF, 2000; Wendel and Malmivuo, 2006. \nbrain and scalp: 0.33 S/m\nskull: 0.022 S/m'))...
                    PropertyGridField('orientation',true,'DisplayName','Orientation free','Description','If true, computes the LF matrix with orientation free dipoles, resulting in a matris Nsensors X 3*Nvertices. If false the LF matrix is computed with dipoles normal to the cortical surface, then the size would be Nsensors X Nvertices')...
                    ];
                hFigure = figure('MenuBar','none','Name','OpenMEEG solver','NumberTitle', 'off','Toolbar', 'none','Units','pixels','Color',obj.container.container.preferences.gui.backgroundColor,...
                    'Resize','off','userData',0);
                position = get(hFigure,'position');
                set(hFigure,'position',[position(1:2) 303 250]);
                hPanel = uipanel(hFigure,'Title','','BackgroundColor','white','Units','pixels','Position',[0 55 303 175],'BorderType','none');
                g = PropertyGrid(hPanel,'Properties', prefObj,'Position', [0 0 1 1]);
                uicontrol(hFigure,'Position',[72 15 70 21],'String','Cancel','ForegroundColor',obj.container.container.preferences.gui.fontColor,...
                    'BackgroundColor',obj.container.container.preferences.gui.buttonColor,'Callback',@cancelCallback);
                uicontrol(hFigure,'Position',[164 15 70 21],'String','Ok','ForegroundColor',obj.container.container.preferences.gui.fontColor,...
                    'BackgroundColor',obj.container.container.preferences.gui.buttonColor,'Callback',@okCallback);
                uiwait(hFigure);
                if ~ishandle(hFigure), return;end
                if ~get(hFigure,'userData'), close(hFigure);return;end
                close(hFigure);
                drawnow;
                val = g.GetPropertyValues();
                conductivity = val.conductivity;
                orientation = val.orientation;
                dispCommand = true;
            end
            
            if dispCommand
                disp('Running:');
                if isa(obj,'coreStreamObject')
                    itemIndex = obj.container.findItem(obj.uuid);
                    fprintf('  mobilab.allStreams.item{%i}.computeLeadFieldBEM( [ %i %i %i ], %i );\n',itemIndex,conductivity(1),conductivity(2),conductivity(3),orientation);
                else fprintf('  obj.computeLeadFieldBEM( [ %i %i %i ], %i );\n',conductivity(1),conductivity(2),conductivity(3),orientation);
                end
            end
            
            if ~exist(obj.surfaces,'file'), error('The file containing the surfaces is missing.');end
            status = system('which om_assemble');
            existOM = ~status;
            if ~existOM
                try
                    mobilab = evalin('base','mobilab');
                    mobilabPath = mobilab.path;
                catch
                    mobilabPath = which('mobilabApplication');
                    if ~isempty(mobilabPath), mobilabPath = fileparts(mobilabPath);
                    else error('OpenMEEG is not intalled. Please download and install the sources you need from https://gforge.inria.fr/frs/?group_id=435.');
                    end
                end
                openmeegDir = [mobilabPath filesep 'dependency' filesep 'openmeeg'];
                
                %---
                % Approach taken from Brainstorm's function bst_openmeeg,
                % Francois Tadel & Alexandre Gramfort, 2011
                %---
                if ~ispc
                    if ismember(computer, {'GLNX86','GLNXA64'}), varname = 'LD_LIBRARY_PATH';
                    else varname = 'DYLD_LIBRARY_PATH';
                    end
                    libpath = getenv(varname);
                    if ~isempty(libpath), libpath = [libpath ':'];end
                    if isempty(strfind(lower(libpath),'openmeeg')), setenv(varname, [libpath openmeegDir]);end
                end
                % Set number of cores used
                try numcores = feature('numcores');
                catch
                    numcores = 4;
                end
                setenv('OMP_NUM_THREADS', num2str(numcores));
                %---
            end
            
            load(obj.surfaces);
            Ns = length(surfData); %#ok
            gTools = geometricTools;            
            rootDir = fileparts(obj.surfaces);
            if isempty(rootDir), rootDir = pwd;end
            binDir = fileparts(which('libmatio.a'));
            [~,rname] = fileparts(tempname);
            headModelGeometry = fullfile(rootDir,[rname '.geom']);
            try %#ok
                copyfile( which('head_model.geom'),headModelGeometry,'f');
                c1 = onCleanup(@()delete(headModelGeometry));
            end
            
            headModelConductivity = fullfile(rootDir,[rname '.cond']);
            fid = fopen(headModelConductivity,'w');
            fprintf(fid,'# Properties Description 1.0 (Conductivities)\n\nAir         0.0\nScalp       %.3f\nBrain       %0.3f\nSkull       %0.3f',...
                conductivity(1),conductivity(3),conductivity(2));
            fclose(fid);
            c2 = onCleanup(@()delete(headModelConductivity));
            
            dipolesFile = fullfile(rootDir,[rname '_dipoles.txt']);
            normalsIn = false;
            [normals,surfData(Ns).faces] = gTools.getSurfaceNormals(surfData(Ns).vertices,surfData(Ns).faces,normalsIn);
            
            normalityConstrained = ~orientation;
            if normalityConstrained, sourceSpace = [surfData(Ns).vertices normals];
            else One = ones(length(normals(:,2)),1);
                Zero = 0*One;
                sourceSpace = [surfData(Ns).vertices One Zero Zero;...
                    surfData(Ns).vertices Zero One Zero;...
                    surfData(Ns).vertices Zero Zero One];
            end
            dlmwrite(dipolesFile, sourceSpace, 'precision', 6,'delimiter',' ')
            c3 = onCleanup(@()delete(dipolesFile));
            
            electrodesFile = fullfile(rootDir,[rname '_elec.txt']);
            dlmwrite(electrodesFile, obj.channelSpace, 'precision', 6,'delimiter',' ')
            c4 = onCleanup(@()delete(electrodesFile));
            
            normalsIn = true;
            brain = fullfile(rootDir,'brain.tri');
            if Ns == 4
                [normals,surfData(3).faces] = gTools.getSurfaceNormals(surfData(3).vertices,surfData(3).faces,normalsIn);
                om_save_tri(brain,surfData(3).vertices,surfData(3).faces,normals)
            else
                [normals,surfData(2).faces] = gTools.getSurfaceNormals(surfData(2).vertices,surfData(2).faces,normalsIn);
                csfSurf = surfData(2);
                csfSurf.vertices     = surfData(2).vertices + 1.05*normals;
                surfData(2).vertices = surfData(2).vertices - 1.05*normals;
                csfSurf.vertices     = gTools.repareIntersectedSurface(surfData(end),csfSurf,1);
                surfData(2).vertices = gTools.repareIntersectedSurface(csfSurf,surfData(2),2);
                surfData(1).vertices = gTools.repareIntersectedSurface(surfData(2),surfData(1),2);
                [normals,csfSurf.faces] = gTools.getSurfaceNormals(csfSurf.vertices,csfSurf.faces,normalsIn);
                om_save_tri(brain,csfSurf.vertices,csfSurf.faces,normals)
            end
            c5 = onCleanup(@()delete(brain));
            
            skull = fullfile(rootDir,'skull.tri');
            [normals,surfData(2).faces] = gTools.getSurfaceNormals(surfData(2).vertices,surfData(2).faces,normalsIn);
            om_save_tri(skull,surfData(2).vertices,surfData(2).faces,normals)
            c6 = onCleanup(@()delete(skull));
            
            head = fullfile(rootDir,'head.tri');
            [normals,surfData(1).faces] = gTools.getSurfaceNormals(surfData(1).vertices,surfData(1).faces,normalsIn);
            om_save_tri(head,surfData(1).vertices,surfData(1).faces,normals)
            c7 = onCleanup(@()delete(head));
            
            hmFile    = fullfile(rootDir,'hm.bin');    c8  = onCleanup(@()delete(hmFile));
            hmInvFile = fullfile(rootDir,'hm_inv.bin');c9  = onCleanup(@()delete(hmInvFile));
            dsmFile   = fullfile(rootDir,'dsm.bin');   c10 = onCleanup(@()delete(dsmFile));
            h2emFile  = fullfile(rootDir,'h2em.bin');  c11 = onCleanup(@()delete(h2emFile));
            lfFile    = fullfile(rootDir,[rname '_LF.mat']);
            
            if ~existOM
                runHere = './';
                wDir = pwd;
                cd(binDir);
            else runHere = '';
            end
            try
                out = system([runHere 'om_assemble -HM "' headModelGeometry '" "' headModelConductivity '" "' hmFile '"']);
                if out, error('An unexpected error occurred running OpenMEEG binaries. Report this to alejandro@sccn.ucsd.edu');end
                
                out = system([runHere 'om_minverser "' hmFile '" "' hmInvFile '"']);
                if out, error('An unexpected error occurred running OpenMEEG binaries. Report this to alejandro@sccn.ucsd.edu');end
                
                out = system([runHere 'om_assemble -DSM "' headModelGeometry '" "' headModelConductivity '" "' dipolesFile '" "' dsmFile '"']);
                if out, error('An unexpected error occurred running OpenMEEG binaries. Report this to alejandro@sccn.ucsd.edu');end
                
                out = system([runHere 'om_assemble -H2EM "' headModelGeometry '" "' headModelConductivity '" "' electrodesFile '" "' h2emFile '"']);
                if out, error('An unexpected error occurred running OpenMEEG binaries. Report this to alejandro@sccn.ucsd.edu');end
                
                out = system([runHere 'om_gain -EEG "' hmInvFile '" "' dsmFile '" "' h2emFile '" "' lfFile '"']);
                if out, error('An unexpected error occurred running OpenMEEG binaries. Report this to alejandro@sccn.ucsd.edu');end
            catch ME
                if strcmp(pwd,binDir), cd(wDir);end
                ME.rethrow;
            end
            if strcmp(pwd,binDir), cd(wDir);end
            if ~exist(lfFile,'file'), error('An unexpected error occurred running OpenMEEG binaries. Report this to alejandro@sccn.ucsd.edu');end
                        
            load(lfFile);
            K = linop;
            clear linop;
            
            %-- Remove extreme values due to numerical instability
            z = zscore(K(:));
            a = 0.00001;
            ind = find(z<norminv(a) | z>norminv(1-a));
            ind_i = setdiff(1:numel(K),ind);
            K(ind) = interp1(ind_i,K(ind_i),ind,'nearest','extrap');
            %--
            
            if exist(lfFile,'file'), delete(lfFile);end
            if exist(obj.leadFieldFile,'file'), delete(obj.leadFieldFile);end
            if isa(obj,'coreStreamObject'), lfFile = fullfile(obj.container.mobiDataDirectory,['lf_' obj.name '_' obj.uuid '_' obj.sessionUUID '.mat']);end
            obj.leadFieldFile = lfFile;
            save(obj.leadFieldFile,'K');
            if isa(obj,'coreStreamObject'), saveProperty(obj,'leadFieldFile',obj.leadFieldFile);end
            disp('Done.')
        end
       %%
        function [sourceSpace,K,L,rmIndices] = getSourceSpace4PEB(obj,structName, rmIndices)
            if isempty(obj.surfaces) || isempty(obj.leadFieldFile), error('Head model or leadfield are missing.');end
            if nargin < 2
                structName = {'Thalamus_L' 'Thalamus_R'};
                disp('Undefined structure to remove. Opening the surface by the Thalamus.')
            end
            if nargin < 3, rmIndices = [];end
            maxNumVertices2rm = 10;
            load(obj.surfaces,'-mat');
            sourceSpace = surfData(end); %#ok
            load(obj.leadFieldFile,'-mat');
            if ~exist('L','var'),
                disp('Computing the Laplacian operator...')
                L = geometricTools.getSurfaceLaplacian(sourceSpace.vertices,sourceSpace.faces);
                save(obj.leadFieldFile,'K','L','-mat')
            end
            
            try 
                [sourceSpace,rmIndices] = obj.removeStructureFromSourceSpace(structName,maxNumVertices2rm, rmIndices);
            catch ME
                warning(ME.message);
                disp('Doing my best to open the surface.')
                n = size(sourceSpace.vertices,1);
                rmIndices = fix(n/2)-maxNumVertices2rm/2:fix(n/2)+maxNumVertices2rm/2;
            end
            dim = size(K); %#ok
            L(rmIndices,:) = [];
            L(:,rmIndices) = [];
            if dim(2)/3 == size(surfData(end).vertices,1) %#ok
                K = reshape(K,[dim(1) dim(2)/3 3]);
                K(:,rmIndices,:) = [];
                % K = permute(K,[1 3 2]);
                K = reshape(K,[dim(1) (dim(2)/3-length(rmIndices))*3]);
                L = kron(eye(3),L);
            else
                K(:,rmIndices) = [];
            end
        end
       %%
        function indices = indices4Structure(obj,structName)
            if nargin < 2, error('Not enough input arguments.');end
            ind = find(ismember(obj.atlas.label,structName));
            if isempty(ind), error('MoBILAB:noStructureMatched','The structure you want to remove is not defined in this atlas.');end
            indices = bsxfun(@eq,obj.atlas.colorTable,ind');
        end
       %%
        function xyz = getCentroidROI(obj,ROInames)
            if nargin < 2, error('Not enough input arguments.');end
            if isempty(obj.atlas) || isempty(obj.surfaces), error('Head model or atlas are missing.');end
            if ~iscell(ROInames), ROInames = {ROInames}; end
            N = length((ROInames));
            xyz = nan(N,3);
            load(obj.surfaces);
            for it=1:N
                try indices = obj.indices4Structure(ROInames{it});
                    xyz(it,:) = mean(surfData(end).vertices(indices,:));
                end
            end
        end
       %%
        function [FP,S] = getForwardProjection(obj,xyz)
            if nargin < 2, error('Not enough input arguments.');end
            if isempty(obj.atlas) || isempty(obj.surfaces), error('Head model or atlas are missing.');end
            if ~exist(obj.surfaces,'file'), error('Head model is missing.');end
            if isempty(obj.leadFieldFile), error('Lead field is missing.');end
            if ~exist(obj.leadFieldFile,'file'), error('Lead field is missing.');end
            load(obj.leadFieldFile);
            if ~exist('K','var'), error('Lead field is missing.');end
            
            load(obj.surfaces);
            [~,~,loc] = geometricTools.nearestNeighbor(surfData(end).vertices,xyz);
            dim = size(K);
            if size(surfData(end).vertices,1) == dim(2)/3, K = reshape(K,[dim(1) dim(2)/3 3]);end
            FP = sum(K(:,loc,:),3);
            S = geometricTools.simulateGaussianSource(surfData(end).vertices,xyz,0.01);
        end
       %%
        function hFigureObj = plotDipoles(obj,xyz,ecd,dipoleLabel,figureTitle)
            if nargin < 2, error('Not enough input arguments.');end
            if isempty(obj.surfaces), error('Head model is missing.');end
            N = size(xyz,1);
            if nargin < 3, ecd = 3*ones(N,3);end
            if isempty(ecd), ecd = 3*ones(N,3);end
            if nargin < 4, dipoleLabel = [];end
            if nargin < 5, figureTitle = '';end
            hFigureObj = equivalentCurrentDipoleViewer(obj,xyz,ecd,dipoleLabel,figureTitle);
        end
        function hFigureObj = plotDipolesForwardProjection(obj,xyz,figureTitle)
            [FP,S] = getForwardProjection(obj,xyz);
            hFigureObj = obj.plotOnModel(S,FP);
        end
       %%
        function [sourceSpace,rmIndices] = removeStructureFromSourceSpace(obj,structName,maxNumVertices2rm, structIndices)
            if isempty(obj.atlas) || isempty(obj.surfaces), error('Head model or atlas are missing.');end
            if nargin < 2, error('Not enough input arguments.');end
            if nargin < 3, maxNumVertices2rm = [];end
            if nargin < 4, structIndices = [];end
            if ~iscell(structName), structName = {structName}; end
            
            load(obj.surfaces,'-mat');
            sourceSpace = surfData(end);%#ok
            
            tmpIndices = indices4Structure(obj,structName);
            if ~any(tmpIndices(:)) && isempty(structIndices),
                error('The structure you want to remove is not defined in this atlas.');
            end
            
            if ~isempty(structIndices)
                % concatenate elements of structIndices into a single column vector
                structIndices = cellfun(@(x)x(:),structIndices,'UniformOutput',false)';
                structIndices = cell2mat(structIndices);
            end
            if ~isempty(maxNumVertices2rm) && any(sum(tmpIndices) > maxNumVertices2rm+1)
                I = [];
                maxNumVertices2rm = fix(maxNumVertices2rm/size(tmpIndices,2));
                for it=1:size(tmpIndices,2)
                    ind = find(tmpIndices(:,it));
                    if length(ind) > maxNumVertices2rm
                        I = [I; ind(1:maxNumVertices2rm)];
                    else
                        I = [I; ind];
                    end
                end
                tmpIndices = I;
            end
            rmIndices = unique_bc([tmpIndices(:) ; structIndices]);
            
            [nVertices,nFaces] = geometricTools.openSurface(sourceSpace.vertices,sourceSpace.faces,rmIndices);
            sourceSpace.vertices = nVertices;
            sourceSpace.faces = nFaces;
        end
       %%
        function saveToFile(obj,file)
            metadata = struct(obj);
            if exist(metadata.surfaces,'file')
                metadata.surfData = load(metadata.surfaces);
            else metadata.surfData = [];
            end
            if exist(metadata.leadFieldFile,'file')
                metadata.leadField = load(metadata.leadFieldFile);
            else metadata.leadField = []; %#ok
            end
            save(file,'metadata','-mat');
        end
        function delete(obj)
            if exist(obj.surfaces,'file')
                [~,filename] = fileparts(obj.surfaces);
                if filename(1) == '.', delete(obj.surfaces);end
            end
            if exist(obj.leadFieldFile,'file')
                [~,filename] = fileparts(obj.leadFieldFile);
                if filename(1) == '.', delete(obj.leadFieldFile);end
            end
        end
    end
    methods(Static)
        function obj = loadFromFile(file)
            metadata = load(file,'-mat');
            if isfield(metadata,'metadata')
                metadata = metadata.metadata;
            end
            if ~isempty(metadata.surfData)
                surfData = metadata.surfData;
                if isfield(surfData,'surfData'), surfData = surfData.surfData;end%#ok
                % [~,filename] = fileparts(tempname);
                % metadata.surfaces = [getHomeDir filesep '.' filename '.mat'];
                metadata.surfaces = [tempname '.mat'];
                save(metadata.surfaces,'surfData');
            end
            if isfield(metadata,'leadField') && ~isempty(metadata.leadField)
                % [~,filename] = fileparts(tempname);
                % metadata.leadFieldFile = [getHomeDir filesep '.' filename '.mat'];
                metadata.leadFieldFile = [tempname '.mat'];
                if isfield(metadata.leadField,'K')
                    K = metadata.leadField.K; %#ok
                else
                    K = metadata.leadField; %#ok
                end
                if isfield(metadata.leadField,'L')
                    L = metadata.leadField.L; %#ok
                    save(metadata.leadFieldFile,'K','L');
                else save(metadata.leadFieldFile,'K');
                end
            else
                metadata.leadFieldFile = [];
            end
            obj = headModel('channelSpace',metadata.channelSpace,'fiducials',metadata.fiducials,'surfaces',metadata.surfaces,...
                'atlas',metadata.atlas,'leadFieldFile',metadata.leadFieldFile,'label',metadata.label);
        end
    end
end

%--
function [elec,labels,fiducials] = readMontage(file)
[eloc, labels] = readlocs(file);
elec = [cell2mat({eloc.X}'), cell2mat({eloc.Y}'), cell2mat({eloc.Z}')];
Nl = length(labels);
count = 1;
lowerLabels = lower(labels);
rmThis = false(Nl,1);
for it=1:Nl
    if ~isempty(strfind(lowerLabels{it},'fidnz')) || ~isempty(strfind(lowerLabels{it},'nasion')) || ~isempty(strfind(lowerLabels{it},'Nz'))
        fiducials.nasion = elec(it,:);
        rmThis(it) = true;
        count = count+1;
    elseif ~isempty(strfind(lowerLabels{it},'fidt9')) || ~isempty(strfind(lowerLabels{it},'lpa')) || ~isempty(strfind(lowerLabels{it},'LPA'))
        fiducials.lpa = elec(it,:);  
        rmThis(it) = true;
        count = count+1;
    elseif ~isempty(strfind(lowerLabels{it},'fidt10')) || ~isempty(strfind(lowerLabels{it},'rpa')) || ~isempty(strfind(lowerLabels{it},'RPA'))
        fiducials.rpa = elec(it,:);
        rmThis(it) = true;
        count = count+1;
    elseif ~isempty(strfind(lowerLabels{it},'fidt10')) || ~isempty(strfind(lowerLabels{it},'vertex'))
        fiducials.vertex = elec(it,:);
        rmThis(it) = true;
        count = count+1;
    end
    if count > 4, break;end
end
elec(rmThis,:) = [];
labels(rmThis) = [];
end


%% unique_bc - unique backward compatible with Matlab versions prior to 2013a
function [C,IA,IB] = unique_bc(A,varargin);

errorFlag = error_bc;

v = version;
indp = find(v == '.');
v = str2num(v(1:indp(2)-1));
if v > 7.19, v = floor(v) + rem(v,1)/10; end;

if nargin > 2
    ind = strmatch('legacy', varargin);
    if ~isempty(ind)
        varargin(ind) = [];
    end;
end;

if v >= 7.14
    [C,IA,IB] = unique(A,varargin{:},'legacy');
    if errorFlag
        [C2,IA2] = unique(A,varargin{:});
        if ~isequal(C, C2) || ~isequal(IA, IA2) || ~isequal(IB, IB2)
            warning('backward compatibility issue with call to unique function');
        end;
    end;
else
    [C,IA,IB] = unique(A,varargin{:});
end
end

%% ismember_bc - ismember backward compatible with Matlab versions prior to 2013a
function [C,IA] = ismember_bc(A,B,varargin);

errorFlag = error_bc;

v = version;
indp = find(v == '.');
v = str2num(v(1:indp(2)-1));
if v > 7.19, v = floor(v) + rem(v,1)/10; end;

if nargin > 2
    ind = strmatch('legacy', varargin);
    if ~isempty(ind)
        varargin(ind) = [];
    end;
end;

if v >= 7.14
    [C,IA] = ismember(A,B,varargin{:},'legacy');
    if errorFlag
        [C2,IA2] = ismember(A,B,varargin{:});
        if (~isequal(C, C2) || ~isequal(IA, IA2))
            warning('backward compatibility issue with call to ismember function');
        end;
    end;
else
    [C,IA] = ismember(A,B,varargin{:});
end
end

%%
function res = error_bc
res = false;
end