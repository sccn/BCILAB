% Defines the class currentSourceViewer for visualization of EEG inverse solutions on a realistic head model. 
% This class is part of MoBILAB software. 
% For more details visit:  https://code.google.com/p/mobilab/
% 
% Author: Alejandro Ojeda, SCCN/INC/UCSD, Jan-2012

classdef equivalentCurrentDipoleViewer < handle
    properties
        hFigure
        hAxes
        streamObj
        hLabels
        hSensors
        hScalp
        hCortex
        hVector
        hDipoles
        dcmHandle
        ecd
        xyz
    end
    methods
        function obj = equivalentCurrentDipoleViewer(streamObj,xyz,ecd,dipoleLabel,figureTitle)
            if nargin < 2, error('Not enough input arguments.');end
            N = size(xyz,1);
            if nargin < 3, ecd = ones(N,3);end
            if nargin < 4, dipoleLabel = [];end
            if nargin < 5, figureTitle = '';end
            if isa(streamObj,'pcdStream'),
                channelLabels = streamObj.parent.label;
            else
                channelLabels = streamObj.channelLabel;
            end
            if size(ecd,2) == 1
                ecd = [ecd,ecd,ecd]/3;
            end
            obj.streamObj = streamObj; 
            load(obj.streamObj.surfaces);
            color = [0.93 0.96 1];
            
            path = fileparts(which('equivalentCurrentDipoleViewer'));
            loc = strfind(path,'mobilab');
            if ~isempty(loc)
                path = path(1:loc+length('mobilab')-1);
            end
            path = fullfile(path,'skin');
            
            try labelsOn  = imread([path filesep 'labelsOn.png']);
                labelsOff = imread([path filesep 'labelsOff.png']);
                sensorsOn = imread([path filesep 'sensorsOn.png']);
                sensorsOff = imread([path filesep 'sensorsOff.png']);
                scalpOn = imread([path filesep 'scalpOn.png']);
                scalpOff = imread([path filesep 'scalpOff.png']);
            catch ME
                ME.rethrow;
            end
            if isa(streamObj,'struct'), visible = 'off';else visible = 'on';end
            obj.hFigure = figure('Menubar','figure','ToolBar','figure','renderer','opengl','Visible',visible,'Color',color,'Name',figureTitle);
            position = get(obj.hFigure,'Position');
            set(obj.hFigure,'Position',[position(1:2) 1.06*position(3:4)]);
            obj.hAxes = axes('Parent',obj.hFigure);         
            toolbarHandle = findall(obj.hFigure,'Type','uitoolbar');
            
            hcb(1) = uitoggletool(toolbarHandle,'CData',labelsOff,'Separator','on','HandleVisibility','off','TooltipString','Labels On/Off','userData',{labelsOn,labelsOff},'State','off');
            set(hcb(1),'OnCallback',@(src,event)rePaint(obj,hcb(1),'labelsOn'),'OffCallback',@(src, event)rePaint(obj,hcb(1),'labelsOff'));
            
            hcb(2) = uitoggletool(toolbarHandle,'CData',sensorsOff,'Separator','off','HandleVisibility','off','TooltipString','Sensors On/Off','userData',{sensorsOn,sensorsOff},'State','off');
            set(hcb(2),'OnCallback',@(src,event)rePaint(obj,hcb(2),'sensorsOn'),'OffCallback',@(src, event)rePaint(obj,hcb(2),'sensorsOff'));
            
            hcb(3) = uitoggletool(toolbarHandle,'CData',scalpOff,'Separator','off','HandleVisibility','off','TooltipString','Scalp On/Off','userData',{scalpOn,scalpOff},'State','off');
            set(hcb(3),'OnCallback',@(src,event)rePaint(obj,hcb(3),'scalpOn'),'OffCallback',@(src, event)rePaint(obj,hcb(3),'scalpOff'));
                        
            obj.dcmHandle = datacursormode(obj.hFigure);
            obj.dcmHandle.SnapToDataVertex = 'off';
            set(obj.dcmHandle,'UpdateFcn',@(src,event)showLabel(obj,event));
            obj.dcmHandle.Enable = 'off';
            hold(obj.hAxes,'on');
            
            obj.hSensors = scatter3(obj.hAxes,obj.streamObj.channelSpace(:,1),obj.streamObj.channelSpace(:,2),...
                obj.streamObj.channelSpace(:,3),'filled','MarkerEdgeColor','k','MarkerFaceColor','y');
            set(obj.hSensors,'Visible','off');
            
            N = length(channelLabels);
            k = 1.1;
            obj.hLabels = zeros(N,1);
            for it=1:N, obj.hLabels(it) = text('Position',k*obj.streamObj.channelSpace(it,:),'String',channelLabels{it},'Parent',obj.hAxes);end
            set(obj.hLabels,'Visible','off');
            
            obj.ecd = ecd;
            obj.xyz = xyz;
            skinColor = [1 0.75 0.65];
            
            % dipoles
            Norm = mean(sqrt(sum(surfData(end).vertices.^2,2)));
            ecd = 0.1*Norm*ecd/norm(ecd);
            for it=1:size(xyz,1)
                [sx,sy,sz] = ellipsoid(xyz(it,1),xyz(it,2),xyz(it,3),ecd(it,1),ecd(it,2),ecd(it,3));
                surf(obj.hAxes,sx,sy,sz,'LineStyle','none','FaceColor','y');
            end
            
            % vectors
            hvx = quiver3(xyz(:,1),xyz(:,2),xyz(:,3),ecd(:,1),0*ecd(:,2),0*ecd(:,3),0.25,'r','LineWidth',2);
            hvy = quiver3(xyz(:,1),xyz(:,2),xyz(:,3),0*ecd(:,1),ecd(:,2),0*ecd(:,3),0.25,'g','LineWidth',2);
            hvz = quiver3(xyz(:,1),xyz(:,2),xyz(:,3),0*ecd(:,1),0*ecd(:,2),ecd(:,3),0.25,'b','LineWidth',2);
            obj.hVector = [hvx;hvy;hvz];
            
            % cortex
            obj.hCortex = patch('vertices',surfData(end).vertices,'faces',surfData(end).faces,'FaceColor',skinColor,...
                'FaceLighting','phong','LineStyle','none','FaceAlpha',0.1,'SpecularColorReflectance',0,...
                'SpecularExponent',50,'SpecularStrength',0.5,'Parent',obj.hAxes);
            camlight(0,180)
            camlight(0,0)
                        
            % scalp
            obj.hScalp = patch('vertices',surfData(1).vertices,'faces',surfData(1).faces,'facecolor',skinColor,...
                'facelighting','phong','LineStyle','none','FaceAlpha',0.5,'Parent',obj.hAxes,'Visible','off');
            
            if ~isempty(dipoleLabel)
                xyz = 1.1*xyz;
                for k=1:length(dipoleLabel), text('Position',xyz(k,:),'String',dipoleLabel{k},'FontSize',14,'Parent',obj.hAxes);end
            end
            view(obj.hAxes,[90 0]);
            hold(obj.hAxes,'off');
            axis(obj.hAxes,'equal','vis3d');
            axis(obj.hAxes,'off')
            if isprop(streamObj,'name'), objName = [streamObj.name ': '];else objName = '';end
            set(obj.hFigure,'Visible',visible,'userData',obj,'Name',['ECD ' objName]);
            rotate3d
            drawnow
        end
        %%
        function rePaint(obj,hObject,opt)
            CData = get(hObject,'userData');
            if isempty(strfind(opt,'Off'))
                set(hObject,'CData',CData{2});
            else
                set(hObject,'CData',CData{1});
            end
            switch opt
                case 'labelsOn',   set(obj.hLabels,'Visible','on');
                case 'labelsOff',  set(obj.hLabels,'Visible','off');
                case 'sensorsOn',  set(obj.hSensors,'Visible','on');
                case 'sensorsOff', set(obj.hSensors,'Visible','off');
                case 'scalpOn',    set(obj.hScalp,'Visible','on');
                case 'scalpOff',   set(obj.hScalp,'Visible','off');
            end
        end
        %%
        function output_txt = showLabel(obj,event_obj)
            persistent DT
            if strcmp(obj.dcmHandle.Enable,'off'),return;end
            if isempty(DT)
                load(obj.streamObj.surfaces);
                vertices = surfData(end).vertices;
                DT = DelaunayTri(vertices(:,1),vertices(:,2),vertices(:,3));
            end
            pos = get(event_obj,'Position');
            loc = nearestNeighbor(DT, pos);
            output_txt = obj.streamObj.atlas.label{obj.streamObj.atlas.colorTable(loc)};
            drawnow
        end
    end
end
