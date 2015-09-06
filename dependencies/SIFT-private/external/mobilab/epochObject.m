classdef epochObject < handle
    properties
        timeStamp
        label
        condition
        eventInterval
        sorting
        subjectID
    end
    properties(Dependent)
        data
    end
    properties(Hidden,SetAccess=protected)
        mmfObj
        binFile
    end
    methods
        function obj = epochObject(varargin)
            if length(varargin) < 6, error('Not enough input arguments.');end
            data      = varargin{1};
            timeStamp = varargin{2}; %#ok
            
            if size(data,1) ~= length(timeStamp) %#ok
                error('The timeStamp vector has to match the first dimension of the matrix data.');
            end
            obj.binFile = tempname;
            fid = fopen(obj.binFile,'w');fwrite(fid,data(:),class(data));fclose(fid);
            obj.mmfObj = memmapfile(obj.binFile,'Format',{class(data) size(data) 'x'},'Writable',true);
            
            obj.timeStamp     = timeStamp(:); %#ok
            obj.label         = varargin{3}(:);
            obj.condition     = varargin{4};
            obj.eventInterval = varargin{5};
            obj.subjectID     = varargin{6};           
            obj.sorting       = 1:size(obj.data,2);
        end
        %%
        function hFigure = eventHistogram(obj,nbins)
            if nargin < 2, nbins = 10;end            
            I = abs(obj.eventInterval) > 10*median(obj.eventInterval);
            strTitle = ['Inter event inteval histogram (' obj.condition ')'];
            hFigure = figure('Name',strTitle,'Color',[0.93 0.96 1]);
            hist(obj.eventInterval(~I),nbins);
            xlabel('Inter event intervals (sec)');
            ylabel('Frequency of occurrence');
            title(strTitle);
        end
        %%
        function data = get.data(obj), data = obj.mmfObj.Data.x;end
        function set.data(obj,data)
            if numel(data) ~= numel(obj.mmfObj.Data.x)
                filename = obj.mmfObj.Filename;
                obj.mmfObj = [];
                fid = fopen(filename,'w');
                fwrite(fid,data(:),class(data(1)));
                fclose(fid);
                obj.mmfObj = memmapfile(filename,'Format',{class(data) size(data) 'x'},'Writable',true);
            else obj.mmfObj.Data.x = data;
            end
        end 
        function delete(obj)
            if exist(obj.binFile,'file'), delete(obj.binFile);end;
        end
        function metadata = saveobj(obj)
            propertyNames = properties(obj);
            for it=1:length(propertyNames), eval(['metadata.' propertyNames{it} ' = obj.' propertyNames{it} ';']);end
            metadata.class = class(obj);
        end
        %%
        function cobj = copyobj(obj) %#ok
            [~,name] = fileparts(tempname);
            tmpfile = fullfile(getHomeDir,['.' name]);
            save(tmpfile,'obj')
            cobj = load(tmpfile,'-mat');
            cobj = cobj.obj;
            delete(tmpfile) 
        end
    end
    methods(Static)
        function obj = loadobj(metadata)
            obj_properties = fieldnames(metadata);
            obj_values = struct2cell(metadata);
            varargIn = cat(1,obj_properties,obj_values);
            Np = length(obj_properties);
            index = [1:Np; Np+1:2*Np];
            varargIn = varargIn(index(:)); %#ok
            obj = eval([metadata.class '(varargIn);']);
        end
    end
    methods(Abstract)
        hFigure = plot(obj)
        rmThis = detectOutliers(obj)
        removeOutliers(obj,rmThis)
        sortingByTrialSimilarity(obj)
    end
end
