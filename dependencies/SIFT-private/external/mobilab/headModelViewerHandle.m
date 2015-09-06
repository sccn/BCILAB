classdef headModelViewerHandle < handle
    properties
        hFigure
        hAxes
        streamObj
        hFiducials
        hLabels
        hSensors
        hScalp
        hSkull
        hCortex
        dcmHandle
        label
        roiLabels = {};
    end
    methods
        function obj = headModelViewerHandle(streamObj,label)
            if nargin < 2
                label = cell(size(streamObj.channelSpace,1),1);
                for it=1:length(label), label{it} = num2str(it);end
            end
            obj.label = label;
            obj.streamObj = streamObj; 
            
            if isa(streamObj,'eeg')
                color = streamObj.container.container.preferences.gui.backgroundColor;
                mobilab = streamObj.container.container;
                path = fullfile(mobilab.path,'skin');
            else
                color = [0.93 0.96 1];
                path = fileparts(which('runmobilab'));
                path = fullfile(path,'skin');
            end
            labelsOn  = imread([path filesep 'labelsOn.png']);
            labelsOff = imread([path filesep 'labelsOff.png']);
            sensorsOn = imread([path filesep 'sensorsOn.png']);
            sensorsOff = imread([path filesep 'sensorsOff.png']);
            scalpOn = imread([path filesep 'scalpOn.png']);
            scalpOff = imread([path filesep 'scalpOff.png']);
            skullOn = imread([path filesep 'skullOn.png']);
            skullOff = imread([path filesep 'skullOff.png']);
            cortexOn = imread([path filesep 'cortexOn.png']);
            cortexOff = imread([path filesep 'cortexOff.png']);
            atlasOn = imread([path filesep 'atlasOn.png']);
            atlasOff = imread([path filesep 'atlasOff.png']);
            roiModeOn  = imread([path filesep 'selectRoiOn.svg.png']);
            roiModeOff  = imread([path filesep 'selectRoiOff.svg.png']);

            obj.hFigure = figure('Menubar','figure','ToolBar','figure','renderer','opengl','Visible','on','Color',color,'Name','Head model');
            position = get(obj.hFigure,'position');
            set(obj.hFigure,'position',[position(1:2) 632 490]);

            obj.hAxes = axes('Parent',obj.hFigure);
            
            toolbarHandle = findall(obj.hFigure,'Type','uitoolbar');
            
            hcb(1) = uitoggletool(toolbarHandle,'CData',labelsOn,'Separator','on','HandleVisibility','off','TooltipString','Labels On/Off','userData',{labelsOn,labelsOff},'State','on');
            set(hcb(1),'OnCallback',@(src,event)rePaint(obj,hcb(1),'labelsOn'),'OffCallback',@(src, event)rePaint(obj,hcb(1),'labelsOff'));
            
            hcb(2) = uitoggletool(toolbarHandle,'CData',sensorsOn,'Separator','off','HandleVisibility','off','TooltipString','Sensors On/Off','userData',{sensorsOn,sensorsOff},'State','on');
            set(hcb(2),'OnCallback',@(src,event)rePaint(obj,hcb(2),'sensorsOn'),'OffCallback',@(src, event)rePaint(obj,hcb(2),'sensorsOff'));
            
            hcb(3) = uitoggletool(toolbarHandle,'CData',scalpOn,'Separator','off','HandleVisibility','off','TooltipString','Scalp On/Off','userData',{scalpOn,scalpOff},'State','on');
            set(hcb(3),'OnCallback',@(src,event)rePaint(obj,hcb(3),'scalpOn'),'OffCallback',@(src, event)rePaint(obj,hcb(3),'scalpOff'));
            
            hcb(4) = uitoggletool(toolbarHandle,'CData',skullOn,'Separator','off','HandleVisibility','off','TooltipString','Skull On/Off','userData',{skullOn,skullOff},'State','on');
            set(hcb(4),'OnCallback',@(src,event)rePaint(obj,hcb(4),'skullOn'),'OffCallback',@(src, event)rePaint(obj,hcb(4),'skullOff'));
            
            hcb(5) = uitoggletool(toolbarHandle,'CData',cortexOn,'Separator','off','HandleVisibility','off','TooltipString','Cortex On/Off','userData',{cortexOn,cortexOff},'State','on');
            set(hcb(5),'OnCallback',@(src,event)rePaint(obj,hcb(5),'cortexOn'),'OffCallback',@(src, event)rePaint(obj,hcb(5),'cortexOff'));
            
            hcb(6) = uitoggletool(toolbarHandle,'CData',atlasOn,'Separator','off','HandleVisibility','off','TooltipString','Atlas On/Off','userData',{atlasOn,atlasOff},'State','on');
            set(hcb(6),'OnCallback',@(src,event)rePaint(obj,hcb(6),'atlasOn'),'OffCallback',@(src, event)rePaint(obj,hcb(6),'atlasOff'));
            
            hcb(7) = uitoggletool(toolbarHandle,'CData',roiModeOn,'Separator','on','HandleVisibility','off','TooltipString','ROI mode On/Off','userData',{roiModeOn,roiModeOff},'State','off');
            set(hcb(7),'OnCallback',@(src,event)rePaint(obj,hcb(7),'roiModeOn'),'OffCallback',@(src, event)rePaint(obj,hcb(7),'roiModeOff'));
            
            obj.dcmHandle = datacursormode(obj.hFigure);
            obj.dcmHandle.SnapToDataVertex = 'off';
            set(obj.dcmHandle,'UpdateFcn',@(src,event)showLabel(obj,event));
            obj.dcmHandle.Enable = 'off';
            
            hold(obj.hAxes,'on');
            
            obj.hSensors = scatter3(obj.hAxes,obj.streamObj.channelSpace(:,1),obj.streamObj.channelSpace(:,2),...
                obj.streamObj.channelSpace(:,3),'filled','MarkerEdgeColor','k','MarkerFaceColor','y');
            N = length(obj.label);
            k = 1.1;
            obj.hLabels = zeros(N,1);
            for it=1:N, obj.hLabels(it) = text('Position',k*obj.streamObj.channelSpace(it,:),'String',obj.label{it},'Parent',obj.hAxes);end
            
            try %#ok
                obj.hFiducials(1) = scatter3(obj.hAxes,obj.streamObj.fiducials.nasion(1),obj.streamObj.fiducials.nasion(2),obj.streamObj.fiducials.nasion(3),'filled',...
                    'MarkerEdgeColor','k','MarkerFaceColor','K');
                obj.hLabels(end+1) = text('Position',1.1*obj.streamObj.fiducials.nasion,'String','Nas','FontSize',12,'FontWeight','bold','Color','k','Parent',obj.hAxes);
                
                obj.hFiducials(2) = scatter3(obj.hAxes,obj.streamObj.fiducials.lpa(1),obj.streamObj.fiducials.lpa(2),obj.streamObj.fiducials.lpa(3),'filled',...
                    'MarkerEdgeColor','k','MarkerFaceColor','K');
                obj.hLabels(end+1) = text('Position',1.1*obj.streamObj.fiducials.lpa,'String','LPA','FontSize',12,'FontWeight','bold','Color','k','Parent',obj.hAxes);
                
                obj.hFiducials(3) = scatter3(obj.hAxes,obj.streamObj.fiducials.rpa(1),obj.streamObj.fiducials.rpa(2),obj.streamObj.fiducials.rpa(3),'filled',...
                    'MarkerEdgeColor','k','MarkerFaceColor','K');
                obj.hLabels(end+1) = text('Position',1.1*obj.streamObj.fiducials.rpa,'String','RPA','FontSize',12,'FontWeight','bold','Color','k','Parent',obj.hAxes);
                
                obj.hFiducials(4) = scatter3(obj.hAxes,obj.streamObj.fiducials.vertex(1),obj.streamObj.fiducials.vertex(2),obj.streamObj.fiducials.vertex(3),'filled',...
                    'MarkerEdgeColor','k','MarkerFaceColor','K');
                obj.hLabels(end+1) = text('Position',1.1*obj.streamObj.fiducials.vertex,'String','Ver','FontSize',12,'FontWeight','bold','Color','k','Parent',obj.hAxes);
                
                obj.hFiducials(5) = scatter3(obj.hAxes,obj.streamObj.fiducials.inion(1),obj.streamObj.fiducials.inion(2),obj.streamObj.fiducials.inion(3),'filled',...
                    'MarkerEdgeColor','k','MarkerFaceColor','K');
                obj.hLabels(end+1) = text('Position',1.1*obj.streamObj.fiducials.inion,'String','Ini','FontSize',12,'FontWeight','bold','Color','k','Parent',obj.hAxes);
            end
            load(obj.streamObj.surfaces);
           
            % cortex
            if ~isempty(obj.streamObj.atlas),
                obj.hCortex = patch('vertices',surfData(end).vertices,'faces',surfData(end).faces,'FaceVertexCData',obj.streamObj.atlas.colorTable,...
                    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
                    'SpecularExponent',50,'SpecularStrength',0.5,'Parent',obj.hAxes);
                camlight(0,180)
                camlight(0,0)
            else
                obj.hCortex = patch('vertices',surfData(end).vertices,'faces',surfData(end).faces,'facecolor','green',...
                    'FaceLighting','gouraud','LineStyle','-','LineWidth',.005,'EdgeColor',[.3 .3 .3],'AmbientStrength',.4,...
                    'FaceAlpha',1,'SpecularColorReflectance',0,'SpecularExponent',50,'SpecularStrength',0.5,'Parent',obj.hAxes);
            end
            % skull
            obj.hSkull = patch('vertices',surfData(2).vertices,'faces',surfData(2).faces,'facecolor','white',...
                'facelighting','phong','LineStyle','none','FaceAlpha',.45,'Parent',obj.hAxes);
            
            % scalp
            skinColor = [1,.75,.65];
            obj.hScalp = patch('vertices',surfData(1).vertices,'faces',surfData(1).faces,'facecolor',skinColor,...
                'facelighting','phong','LineStyle','none','FaceAlpha',.25,'Parent',obj.hAxes);
            
            camlight(0,180)%,'infinite') gouraud
            camlight(0,0)
            view(obj.hAxes,[90 0]);
            hold(obj.hAxes,'off');
            axis(obj.hAxes,'equal','vis3d');
            axis(obj.hAxes,'off')
            set(obj.hFigure,'Visible','on','userData',obj);
            rotate3d
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
                case 'labelsOn'
                    set(obj.hLabels,'Visible','on');
                case 'labelsOff'
                    set(obj.hLabels,'Visible','off');
                case 'sensorsOn'
                    set(obj.hSensors,'Visible','on');
                    set(obj.hFiducials,'Visible','on');
                case 'sensorsOff'
                    set(obj.hSensors,'Visible','off');
                    set(obj.hFiducials,'Visible','off');
                case 'scalpOn'
                    set(obj.hScalp,'Visible','on');
                case 'scalpOff'
                    set(obj.hScalp,'Visible','off');
                case 'skullOn'
                    set(obj.hSkull,'Visible','on');
                case 'skullOff'
                    set(obj.hSkull,'Visible','off');
                case 'cortexOn'
                    set(obj.hCortex,'Visible','on');
                case 'cortexOff'
                    set(obj.hCortex,'Visible','off');
                case 'atlasOn'
                    set(obj.hCortex,'FaceVertexCData',obj.streamObj.atlas.colorTable,'LineStyle','none','FaceColor','interp');
                case 'atlasOff'
                    set(obj.hCortex,'FaceColor',[ 0 1 0],'LineStyle','-');
                case 'roiModeOn'
                    set(obj.dcmHandle,'UpdateFcn',@(src,event)showLabel(obj,event,true));
                    obj.roiLabels = {};
                case 'roiModeOff'
                    set(obj.dcmHandle,'UpdateFcn',@(src,event)showLabel(obj,event,false));
            end
        end
        %%
        function output_txt = showLabel(obj,event_obj, storeLabel)
            persistent DT
            if strcmp(obj.dcmHandle.Enable,'off'),return;end
            if nargin < 3, storeLabel = true;end
            if isempty(DT)
                load(obj.streamObj.surfaces);
                vertices = surfData(end).vertices;
                %DT = DelaunayTri(vertices(:,1),vertices(:,2),vertices(:,3));
                DT = delaunayTriangulation(vertices(:,1),vertices(:,2),vertices(:,3));
            end
            pos = get(event_obj,'Position');
            loc = nearestNeighbor(DT, pos);
            output_txt = obj.streamObj.atlas.label{obj.streamObj.atlas.colorTable(loc)};
            
            % store label in cell array
            if storeLabel, obj.roiLabels = unique([obj.roiLabels(:); output_txt]);end
        end
    end
end