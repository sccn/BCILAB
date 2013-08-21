function vis_filtered(varargin)
% Display a filtered representation of a BCILAB stream.
%
% Keyboard shortcuts:
%   [up arrow]   : increase the y scale of the time series
%   [down arrow] : decrease the y scale of the time series
%   [right arrow]: increase the displayed time range
%   [left arrow] : decrease the displayed time range
%   [page up]    : go up by one page of channels
%   [page down]  : go down by one page of channels
%
% In:
%   StreamName : Stream to display. The name of the stream that you would like to display.
%
%   TimeScale : Initial time scale in seconds. The time range of the display window;
%               can be changed with keyboard shortcuts (see help). Default=5
%
%   DataScale : Initial scale of the data. The scale of the data, in units between horizontal lines;
%               can be changed with keyboard shortcuts (see help). Default=150
%
%   ChannelRange : Channels to display. The channel range to display. Default=[1:32]
%
%   SamplingRate : Sampling rate for display. This is the sampling rate that is used for plotting, in Hz;
%                  for faster drawing. Default=100
%
%   RefreshRate : Refresh rate for display. This is the rate at which the graphics are updated, in Hz.
%                 Default=10
%
%   Rereference : Apply common-average re-referencing to the data. Useful for noisy EEG streams.
%                 Default=false
%
%   PageOffset : Channel page offset. Allows to flip forward or backward pagewise through the displayed channels.
%                Default=0
%
%   Position : Figure position. Allows to script the position at which the figures should appear.
%              This is a 4-element vector of the form [X-offset,Y-offset,Width,Height]
%              with all values in pixes.
%              Default=[]
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-07-10
%
%                                uses portions of vis_dataStreamViewer
%                                (c) 2012 by Tim Mullen


streamnames = {};
pipelinenames = {};

% list the streams and pipelines in the workspace
vars = evalin('base','whos');
for vname = {vars.name}
    vname = vname{1};
    var = evalin('base',vname);
    if isfield(var,'tracking') && isfield(var.tracking,'online_expression')
        pipelinenames{end+1} = vname; end
    if all(isfield(var,{'smax','buffer'}))
        streamnames{end+1} = vname; end
end
if isempty(pipelinenames)
    error('No filter pipelines are defined in the workspace. Tip: any preprocessed data set can be used as a reference definition of a filter pipeline'); end
if isempty(streamnames)
    error('No streams are currently running. Tip: you can start a new stream by executing one of the run_read* functions, e.g. run_readlsl.'); end

% handle input arguments
opts = arg_define(varargin, ...
    arg({'streamname','StreamName'},streamnames{1},streamnames,'BCILAB stream to process. The name of the stream that you would like to display.'), ...
    arg({'pipelinename','PipelineName'},pipelinenames{1},pipelinenames,'Filter pipeline to apply. The name of the filter pipeline that you would like to apply.'), ...
    arg({'datafield','DataField'},'data',{'data','icaact','srcpot'},'Data field to plot'), ...
    arg_nogui({'visname','VisualizationName'},[],[],'Name of the visualization. This determines the workspace variable that holds the pipeline and settings.'), ...
    arg({'bufferrange','BufferRange'},10,[],'Maximum time range to buffer. Imposes an upper limit on what can be displayed.'), ...
    arg({'timerange','TimeRange'},5,[],'Initial time range in seconds. The time range of the display window; can be changed with keyboard shortcuts (see help).'), ...
    arg({'datascale','DataScale'},150,[],'Initial scale of the data. The scale of the data, in units between horizontal lines; can be changed with keyboard shortcuts (see help).'), ...
    arg({'channelrange','ChannelRange'},1:32,[],'Channels to display. The channel range to display.'), ...
    arg({'samplingrate','SamplingRate'},100,[],'Sampling rate for display. This is the sampling rate that is used for plotting; for faster drawing.'), ...
    arg({'refreshrate','RefreshRate'},5,[],'Refresh rate for display. This is the rate at which the graphics are updated.'), ...
    arg({'reref','Rereference'},false,[],'Common average reference. Enable this to view the data with a common average reference filter applied.'), ...
    arg_nogui({'pageoffset','PageOffset'},0,[],'Channel page offset. Allows to flip forward or backward pagewise through the displayed channels.'), ...
    arg_nogui({'position','Position'},[],[],'Figure position. Allows to script the position at which the figures should appear.'));


if isempty(varargin)
    % bring up GUI dialog if no arguments were passed
    arg_guidialog;
else
    % fix up some arguments
    opts.bufferrange = max(opts.bufferrange,opts.timerange)*1.1;    

    % init shared handles
    [fig,ax,lines,axscale,scaleLine] = deal([]);
    
    % create a unique name for this visualization
    if ~isempty(opts.visname)
        visname = opts.visname;
    else
        taken = evalin('base','whos(''vis_*'')');
        visname = genvarname(['vis_' opts.streamname],{taken.name});
    end
        
    % and store the data in a workspace variable
    srate = evalin('base',[opts.streamname '.srate']);
    visinfo.pipeline = onl_newpipeline(evalin('base',opts.pipelinename),opts.streamname);
    visinfo.opts = opts;
    assignin('base',visname,visinfo);
    
    % create the figure
    create_figure(opts);
    
    % set up a timer that updates the visualization
    th = timer('Period', 1.0/opts.refreshrate,'ExecutionMode','fixedRate','TimerFcn',@on_timer,'StartDelay',0.2,'Tag',[visname '_timer']);
    start(th);
end


    % === nested functions (sharing some handles with each other) ===

    % create a new figure and axes
    function create_figure(opts)
        fig = figure('Tag',['Fig' visname],'Name',['Stream:' opts.streamname ';Filter:' opts.pipelinename],'CloseRequestFcn','delete(gcbf)', ...
            'KeyPressFcn',@(varargin)on_key(varargin{2}.Key));
        if ~isempty(opts.position)
            set(fig,'Position',opts.position,'Units','Normalized'); end
        ax = axes('Parent',fig, 'Tag','LSLViewer', 'YDir','reverse');
    end

    function on_timer(varargin)
        try 
            % check if the buffer is still there
            if evalin('base',['exist(''' opts.streamname ''',''var'') && exist(''' visname ''',''var'')'])
                
                
                % === update buffer contents (happens in the base workspace) ===
                
                % pull a new chunk from the stream and process it
                visinfo = evalin('base',visname);
                [stream,visinfo.pipeline] = onl_filtered(visinfo.pipeline,1.1*visinfo.opts.timerange*srate);
                assignin('base',visname,visinfo);
                
                % === data post-processing for plotting ===
                
                % determine channels and samples to display
                plotchans = visinfo.opts.channelrange + visinfo.opts.pageoffset*length(visinfo.opts.channelrange);
                if isempty(plotchans)
                    plotchans = 1:size(stream.(opts.datafield),1);
                else
                    plotchans = intersect(1:size(stream.(opts.datafield),1),plotchans);
                end
                plotdata = stream.(opts.datafield)(plotchans, round(1 : stream.srate/visinfo.opts.samplingrate : end));
                plottime = linspace(stream.xmin,stream.xmax,size(plotdata,2));
                
                % re-reference
                if visinfo.opts.reref
                    plotdata = bsxfun(@minus,plotdata,mean(plotdata)); end
                
                % zero-mean
                plotdata = bsxfun(@minus, plotdata, mean(plotdata,2));
                
                % arrange for plotting
                plotoffsets = (0:size(plotdata,1)-1)*visinfo.opts.datascale;
                plotdata = bsxfun(@plus, plotdata', plotoffsets);
                
                
                % === actual drawing ===
                
                % draw the block contents...
                if ~isempty(plotdata)
                    if ~exist('lines','var') || isempty(lines) || length(lines) ~= size(plotdata,2)
                        lines = plot(ax,plottime,plotdata);
                        title(ax,visinfo.opts.streamname);
                        xlabel(ax,'Time (sec)','FontSize',12);
                        ylabel(ax,'Activation','FontSize',12);
                        % update the axis tickmarks
                    else
                        for k=1:length(lines)
                            set(lines(k),'Ydata',plotdata(:,k));
                            set(lines(k),'Xdata',plottime);
                        end
                    end
                
                end
                
                % update the data scale
                set(ax, 'YTick',plotoffsets, 'YTickLabel',{stream.chanlocs(plotchans).labels});
                axis(ax,[-visinfo.opts.timerange 0 -visinfo.opts.datascale size(plotdata,2)*visinfo.opts.datascale + visinfo.opts.datascale]);                

                drawnow;
            else
                try 
                    disp(['Deleting timer ' get(th,'Tag') '.']);
                catch e
                    disp('Deleting timer.');
                end
                % delete the timer
                warning off MATLAB:timer:deleterunning
                delete(th);
            end
        catch e
            if isempty(findobj('Tag',['Fig' visname]))
                disp('Figure was closed.');
            else
                disp('An error occurred during the stream viewer update: ');
                hlp_handleerror(e);
            end
            warning off MATLAB:timer:deleterunning
            delete(th);
        end
    end

    function on_key(key)
        visinfo = evalin('base',visname);
        switch lower(key)
            case 'uparrow'
                % decrease datascale
                visinfo.opts.datascale = visinfo.opts.datascale*0.9;
            case 'downarrow'
                % increase datascale
                visinfo.opts.datascale = visinfo.opts.datascale*1.1;
            case 'rightarrow'
                % increase timerange
                visinfo.opts.timerange = visinfo.opts.timerange*1.1;                
            case 'leftarrow'
                % decrease timerange
                visinfo.opts.timerange = visinfo.opts.timerange*0.9;                
            case 'pagedown'
                % shift display page offset down
                visinfo.opts.pageoffset = visinfo.opts.pageoffset+1;                
            case 'pageup'
                % shift display page offset up
                visinfo.opts.pageoffset = visinfo.opts.pageoffset-1;
        end
        assignin('base',visname,visinfo);
    end
    
end

