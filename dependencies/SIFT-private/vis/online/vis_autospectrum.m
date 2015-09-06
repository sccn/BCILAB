function [fig opts] = vis_autospectrum(varargin)

% 
% hpower = figure('DeleteFcn',@(varargin)evalin('base','opts.plotSpectrum=false;'));
% haxpwr = axes;
% pwr = zeros(length(PowerChannelsIdx),length(Frequencies));
% 
% 
% 
% for ch=1:length(PowerChannelsIdx)
%     pwr(ch,:) = squeeze(EEG.CAT.Conn.S(PowerChannelsIdx(ch),PowerChannelsIdx(ch),:,:));
% end
% offset = 10*cumsum(ones(length(PowerChannelsIdx),1));
% h=plot(haxpwr,bsxfun(@plus,pwr,offset)');
% set(h,'xdata',EEG.CAT.Conn.freqs);
% axis(haxpwr,'tight');
% 


% Display spectrum
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

% handle input arguments
opts = arg_define(varargin, ...
    arg_norep({'stream','EEG'},[],[],'EEGLAB dataset. Must contain CAT.Conn. Can also the name of a variable in the workspace'), ...
    arg({'datascale','DataScale'},150,[],'Initial scale of the data. The scale of the data, in units between horizontal lines; can be changed with keyboard shortcuts (see help).'), ...
    arg({'freqrange','FrequencyRange'},[],[],'Frequency range. Default is all available freqs','shape','row'), ...
    arg({'channelrange','ChannelRange'},[],[],'Channels to display. The channel range to display. Default is all available channels','shape','row'), ...
    arg({'pageoffset','PageOffset'},0,[],'Channel page offset. Allows to flip forward or backward pagewise through the displayed channels.'), ...
    arg({'position','Position'},[],[],'Figure position. Allows to script the position at which the figures should appear.'), ...
    arg_nogui({'figHandle','FigureHandle'},[],[],'Handle to existing figure to update. Otherwise new figure is created'));

fig = [];

if isempty(varargin)
    % bring up GUI dialog if no arguments were passed
    arg_guidialog;
end

if isempty(opts.stream)
    return;
end

if ischar(opts.stream)
    % grab stream from base
    try
        opts.stream = evalin('base',opts.stream);
    catch err
        hlp_handleerror(err);
        return;
    end
end
    
if isempty(opts.figHandle) || ~ishandle(opts.figHandle)
    % init shared handles
    [fig,ax,lines] = deal([]);

    % create the figure
    create_figure(opts);
else
    fig = opts.figHandle;
    ax = findobj(fig,'type','axes');
    lines = findobj(ax,'tag','lines');
end
    

stream = opts.stream;

% make sure dataset is consistent
res = hlp_checkeegset(stream,{'conn'});
if ~isempty(res)
    error(res);
end

if ~isfield(stream.CAT.Conn,'S')
    fprintf('vis_autospectrum: no spectrum data available\n');
    return;
end

if isempty(opts.channelrange)
    opts.channelrange = 1:size(stream.CAT.Conn.S,1);
end


% === nested functions (sharing some handles with each other) ===

% create a new figure and axes
function create_figure(opts)
    fig = figure('Tag','Fig_SpecViewer','Name','Spectrum Viewer','CloseRequestFcn','delete(gcbf)', ...
        'KeyPressFcn',@(varargin)on_key(varargin{2}.Key));
    if ~isempty(opts.position)
        set(fig,'Position',opts.position); end
    ax = axes('Parent',fig, 'Tag','SpecViewer', 'YDir','reverse');
end



try
    % === data post-processing for plotting ===
    
    % determine channels and samples to display
    plotchans = opts.channelrange + opts.pageoffset*length(opts.channelrange);
    if isempty(plotchans)
        plotchans = 1:stream.CAT.nbchan;
    else
        plotchans = intersect_bc(1:stream.CAT.nbchan,plotchans);
    end
    
    % get range of frequencies
    [freqidx]   = getindex(stream.CAT.Conn.freqs,opts.freqrange);
    plotfreqidx = freqidx(1):freqidx(2);
    
    % extract the autospectra
    % TODO: elim. for loop by indexing directly into a reshaped Conn.S
    % sz = size(tmp.CAT.Conn.S);
    % foo = reshape(stream.CAT.Conn.S,[sz(1)^2 sz(3) sz(4)]);
    % i=plotchans; isequal(foo((i-1).*sz(1)+i,:,:),squeeze(stream.CAT.Conn.S(i,i,:,:)))
    plotpwr = zeros(length(plotchans),length(plotfreqidx));
    for ch=1:length(plotchans)
        plotpwr(ch,:) = squeeze(stream.CAT.Conn.S(plotchans(ch),plotchans(ch),plotfreqidx,:));
    end
    
    plotfreq = stream.CAT.Conn.freqs(plotfreqidx);
    
    % zero-mean
    plotpwr = bsxfun(@minus, plotpwr, mean(plotpwr,2));
    
    % arrange for plotting
    plotoffsets = (0:size(plotpwr,1)-1)*opts.datascale;
    plotpwr = bsxfun(@plus, plotpwr', plotoffsets);
    
    
    % === actual drawing ===
    
    % draw the block contents...
    if ~isempty(plotpwr) 
        if ~exist('lines','var') || isempty(lines) || length(lines) ~= size(plotpwr,1);
            lines = plot(ax,plotfreq,plotpwr,'tag','lines');
            title(ax,'VAR Autospectra');
            xlabel(ax,'Frequency (Hz)','FontSize',12);
            ylabel(ax,'Power','FontSize',12);
        else
            for k=1:length(lines)
                set(lines(k),'Ydata',plotpwr(:,k));
                set(lines(k),'Xdata',plotfreq);
            end
        end
        
        % update the axis limit and tickmarks
        axis(ax,[plotfreq(1) plotfreq(end) -opts.datascale size(plotpwr,2)*opts.datascale + opts.datascale]);
        set(ax, 'YTick',plotoffsets, 'YTickLabel',stream.CAT.curComponentNames(plotchans));
    end
    
    drawnow;
    
catch e
    if isempty(findobj('Tag','Fig_SpecViewer'))
        disp('Figure was closed.');
    else
        disp('An error occurred during the SpecViewer update: ');
        hlp_handleerror(e);
    end
end


function on_key(key)
specOpts = evalin('base','specOpts');
switch lower(key)
    case 'uparrow'
        % decrease datascale
        specOpts.datascale = specOpts.datascale*0.9;
    case 'downarrow'
        % increase datascale
        specOpts.datascale = specOpts.datascale*1.1;
    case 'rightarrow'
        % increase timerange
%         specOpts.timerange = specOpts.timerange*1.1;
    case 'leftarrow'
        % decrease timerange
%         specOpts.timerange = specOpts.timerange*0.9;
    case 'pagedown'
        % shift display page offset down
        specOpts.pageoffset = specOpts.pageoffset+1;
    case 'pageup'
        % shift display page offset up
        specOpts.pageoffset = specOpts.pageoffset-1;
end
assignin('base','specOpts',specOpts);
end

end

