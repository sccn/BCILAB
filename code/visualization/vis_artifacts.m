function [h_old,h_new] = vis_artifacts(new,old,varargin)
% vis_artifacts(NewEEG,OldEEG,Options...)
% Display the artifact rejections done by any of the artifact cleaning functions.
%
% Keyboard Shortcuts:
%   [n] : display just the new time series
%   [o] : display just the old time series
%   [b] : display both time series super-imposed
%   [d] : display the difference between both time series
%   [+] : increase signal scale
%   [-] : decrease signal scale
%   [*] : expand time range
%   [/] : reduce time range
%   [h] : show/hide slider
%   [e] : toggle events
%   [l] : toggle event legend
%
% In:
%   NewEEG     : cleaned continuous EEG data set
%   OldEEG     : original continuous EEG data set
%   Options... : name-value pairs specifying the options, with names:
%                'YRange' : y range of the figure that is occupied by the signal plot
%                'YScaling' : distance of the channel time series from each other in std. deviations
%                'WindowLength : window length to display
%                'NewColor' : color of the new (i.e., cleaned) data
%                'OldColor' : color of the old (i.e., uncleaned) data
%                'HighpassOldData' : whether to high-pass the old data if not already done
%                'ScaleBy' : the data set according to which the display should be scaled, can be 
%                            'old' or 'new' (default: 'new')
%                'ChannelSubset' : optionally a channel subset to display
%                'TimeSubet' : optionally a time subrange to display
%                'DisplayMode' : what should be displayed: 'both', 'new', 'old', 'diff'
%                'ShowSetname' : whether to display the dataset name in the title
%                'EqualizeChannelScaling' : optionally equalize the channel scaling
%                See also code for more options.
%
% Notes:
%   This function is primarily meant for testing purposes and is not a particularly high-quality
%   implementation.
%
% Examples:
%  vis_artifacts(clean,raw)
%
%  % display only a subset of channels
%  vis_artifacts(clean,raw,'ChannelSubset',1:4:raw.nbchan);
%
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-07-06
%
%                                relies on the findjobj() function by Yair M. Altman.

% Copyright (C) Christian Kothe, SCCN, 2012, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

have_signallegend = false;
have_eventlegend = false;

if nargin < 2
    old = new; 
elseif ischar(old)
    varargin = [{old} varargin];
    old = new;
end

% parse options
opts = hlp_varargin2struct(varargin, ...
    {'yrange','YRange'}, [0.05 0.95], ...       % y range of the figure occupied by the signal plot
    {'yscaling','YScaling'}, 3.5, ...           % distance of the channel time series from each other in std. deviations
    {'wndlen','WindowLength'}, 10, ...          % window length to display
    {'newcol','NewColor'}, [0 0 0.5], ...       % color of the new (i.e., cleaned) data
    {'oldcol','OldColor'}, [1 0 0], ...         % color of the old (i.e., uncleaned) data
    {'highpass_old','HighpassOldData'},true, ...% whether to high-pass the old data if not already done
    {'show_removed_portions','ShowRemovedPortions'},true, ...% whether to show removed data portions (if only one set is passed in)
    {'show_events','ShowEvents'},true, ...      % whether to show events
    {'show_eventlegend','ShowEventLegend'},false, ...  % whether to show a legend for the currently visible events
    {'scale_by','ScaleBy'},'allnew',...         % the data set according to which the display should be scaled (can be allold, allnew, wndold, or wndnew)
    {'channel_subset','ChannelSubset'},[], ...  % optionally a channel subset to display
    {'time_subset','TimeSubset'},[],...         % optionally a time subrange to display
    {'display_mode','DisplayMode'},'both',...   % what should be displayed: 'both', 'new', 'old', 'diff'
    {'show_setname','ShowSetname'},true,...     % whether to display the dataset name in the title
    {'line_spec','LineSpec'},'-',...            % line style for plotting
    {'line_width','LineWidth'},0.5,...          % line width
    {'add_legend','AddLegend'},false,...        % add a signal legend
    {'equalize_channel_scaling','EqualizeChannelScaling'},false);  % optionally equalize the channel scaling

% ensure that the data are not epoched and expand the rejections with NaN's (now both should have the same size)
if opts.show_removed_portions
    new = expand_rejections(to_continuous(new));
    old = expand_rejections(to_continuous(old));
end
new.chanlocs = old.chanlocs;

% correct for filter delay
if isfield(new.etc,'filter_delay')
    new.data = new.data(:,[1+round(new.etc.filter_delay*new.srate):end end:-1:(end+1-round(new.etc.filter_delay*new.srate))]); end
if isfield(old.etc,'filter_delay')
    old.data = old.data(:,[1+round(old.etc.filter_delay*old.srate):end end:-1:(end+1-round(old.etc.filter_delay*old.srate))]); end

% make sure that the old data is high-passed the same way as the new data
if opts.highpass_old && isfield(new.etc,'clean_drifts_kernel') && ~isfield(old.etc,'clean_drifts_kernel')
    old.data = old.data';
    for c=1:old.nbchan
        old.data(:,c) = filtfilt_fast(new.etc.clean_drifts_kernel,1,old.data(:,c)); end
    old.data = old.data';
end

if isscalar(opts.line_width)
    opts.line_width = [opts.line_width opts.line_width]; end

% optionally pick a subrange to work on
if ~isempty(opts.channel_subset)
    old = pop_select(old,'channel',opts.channel_subset);
    new = pop_select(new,'channel',opts.channel_subset);
end

if ~isempty(opts.time_subset)
    old = pop_select(old,'time',opts.time_subset);
    new = pop_select(new,'time',opts.time_subset);
end

if opts.equalize_channel_scaling    
    rescale = 1./mad(old.data,[],2);
    new.data = bsxfun(@times,new.data,rescale);
    old.data = bsxfun(@times,old.data,rescale);
end

% generate event colormap
if ~isempty(old.event)
    opts.event_colormap = gen_colormap(old.event,'jet'); end

% calculate whole-data scale
old_iqr = 2*mad(old.data',1)';
old_iqr(isnan(old_iqr)) = deal(mean(old_iqr(~isnan(old_iqr))));
new_iqr = 2*mad(new.data',1)';
new_iqr(isnan(new_iqr)) = deal(mean(new_iqr(~isnan(new_iqr))));

% create figure & slider
lastPos = 0;
hFig = figure('ResizeFcn',@on_window_resized,'KeyPressFcn',@(varargin)on_key(varargin{2}.Key)); hold; axis();
hAxis = gca;
hSlider = uicontrol('style','slider','KeyPressFcn',@(varargin)on_key(varargin{2}.Key)); on_resize();
jSlider = findjobj(hSlider);
jSlider.AdjustmentValueChangedCallback = @on_update;

% do the initial update
on_update();


    function repaint(relPos,moved)
        % repaint the current data
        
        % if this happens, we are maxing out MATLAB's graphics pipeline: let it catch up
        if relPos == lastPos && moved
            return; end
        
        % axes
        cla;
        
        % compute pixel range from axis properties
        xl = get(hAxis,'XLim');
        yl = get(hAxis,'YLim');
        fp = get(hFig,'Position');
        ap = get(hAxis,'Position');
        pixels = (fp(3))*(ap(3)-ap(1));
        ylr = yl(1) + opts.yrange*(yl(2)-yl(1));
        channel_y = (ylr(2):(ylr(1)-ylr(2))/(size(new.data,1)-1):ylr(1))';
        
        % compute sample range
        wndsamples = opts.wndlen * new.srate;
        pos = floor((size(new.data,2)-wndsamples)*relPos);
        wndindices = 1 + floor(0:wndsamples/pixels:(wndsamples-1));
        wndrange = pos+wndindices;
        
        oldwnd = old.data(:,wndrange);
        newwnd = new.data(:,wndrange);
        switch opts.scale_by
            case 'allnew'                
                iqrange = new_iqr;
            case 'allold'
                iqrange = old_iqr;
            case {'wndnew','new'}
                iqrange = mad(newwnd',1)';
                iqrange(isnan(iqrange)) = mad(oldwnd(isnan(iqrange),:)',1)';
            case {'wndold','old'}
                iqrange = mad(oldwnd',1)';
            otherwise
                error('Unsupported scale_by option.');
        end
        scale = ((ylr(2)-ylr(1))/size(new.data,1)) ./ (opts.yscaling*iqrange); scale(~isfinite(scale)) = 0;
        scale(scale>median(scale)*3) = median(scale);
        scale = repmat(scale,1,length(wndindices));
                
        % draw
        if opts.show_setname
            tit = sprintf('%s - ',[old.filepath filesep old.filename]);
        else
            tit = '';
        end
        
        if ~isempty(wndrange)
            tit = [tit sprintf('[%.1f - %.1f]',new.xmin + (wndrange(1)-1)/new.srate, new.xmin + (wndrange(end)-1)/new.srate)];        
        end
        
        xrange = xl(1):(xl(2)-xl(1))/(length(wndindices)-1):xl(2);
        yoffset = repmat(channel_y,1,length(wndindices));
        switch opts.display_mode            
            case 'both'                
                title([tit '; superposition'],'Interpreter','none');
                h_old = plot(xrange, (yoffset + scale.*oldwnd)','Color',opts.oldcol,'LineWidth',opts.line_width(1));
                h_new = plot(xrange, (yoffset + scale.*newwnd)','Color',opts.newcol,'LineWidth',opts.line_width(2));
            case 'new'
                title([tit '; cleaned'],'Interpreter','none');
                plot(xrange, (yoffset + scale.*newwnd)','Color',opts.newcol,'LineWidth',opts.line_width(2));
            case 'old'
                title([tit '; original'],'Interpreter','none');
                plot(xrange, (yoffset + scale.*oldwnd)','Color',opts.oldcol,'LineWidth',opts.line_width(1));
            case 'diff'
                title([tit '; difference'],'Interpreter','none');
                plot(xrange, (yoffset + scale.*(oldwnd-newwnd))','Color',opts.newcol,'LineWidth',opts.line_width(1));
        end
        
        % also plot events
        if opts.show_events && ~isempty(old.event)
            evtlats = [old.event.latency];
            evtindices = find(evtlats>wndrange(1) & evtlats<wndrange(end));
            if ~isempty(evtindices)
                evtpos = xl(1) + (evtlats(evtindices)-wndrange(1))/wndsamples*(xl(2)-xl(1));                
                evttypes = {old.event(evtindices).type};
                % for each visible type
                visible_types = unique(evttypes);
                handles = [];
                labels = {};
                for ty = visible_types(:)'
                    % plot line instances in the right color
                    curtype = ty{1};
                    curcolor = opts.event_colormap.values(strcmp(opts.event_colormap.keys,curtype),:);
                    matchpos = strcmp(evttypes,curtype);
                    h = line([evtpos(matchpos);evtpos(matchpos)],repmat([0;1],1,nnz(matchpos)),'Color',curcolor);
                    handles(end+1) = h(1);
                    labels{end+1} = curtype;
                end
                if opts.show_eventlegend
                    legend(handles,labels,'Location','NorthWest'); 
                    have_eventlegend = true;
                elseif have_eventlegend
                    legend off;
                    have_eventlegend = false;
                end
            end
        end        
        axis([0 1 0 1]);
        
        if opts.add_legend && ~have_signallegend
            legend([h_old(1);h_new(1)],'Original','Corrected');
            have_signallegend = 1;
        end
        drawnow;


        lastPos = relPos;
    end


    function on_update(varargin)
        % slider moved
        repaint(get(hSlider,'Value'),~isempty(varargin));
    end

    function on_resize(varargin)
        % adapt/set the slider's size
        wPos = get(hFig,'Position');
        if ~isempty(hSlider)
            try
                set(hSlider,'Position',[20,20,wPos(3)-40,20]);
            catch,end
            on_update;
        end
    end

    function on_window_resized(varargin)
        % window resized
        on_resize();
    end

    function EEG = to_continuous(EEG)
        % convert an EEG set to continuous if currently epoched
        if ndims(EEG.data) == 3
            EEG.data = EEG.data(:,:);
            [EEG.nbchan,EEG.pnts,EEG.trials] = size(EEG.data);
        end
    end

    function EEG = expand_rejections(EEG)
        % reformat the new data so that it can be super-imposed with the old data
        [EEG.nbchan,EEG.pnts] = size(EEG.data);
        if ~isfield(EEG.etc,'clean_channel_mask')
            EEG.etc.clean_channel_mask = true(1,EEG.nbchan); end
        if ~isfield(EEG.etc,'clean_sample_mask')
            EEG.etc.clean_sample_mask = true(1,EEG.pnts); end
        tmpdata = nan(length(EEG.etc.clean_channel_mask),length(EEG.etc.clean_sample_mask));
        tmpdata(EEG.etc.clean_channel_mask,EEG.etc.clean_sample_mask) = EEG.data;
        EEG.data = tmpdata;
        [EEG.nbchan,EEG.pnts] = size(EEG.data);
    end

    function on_key(key)
        switch lower(key)
            case {'add','+'}
                % decrease datascale
                opts.yscaling = opts.yscaling*0.9;
            case {'subtract','-'}
                % increase datascale
                opts.yscaling = opts.yscaling*1.1;
            case {'multiply','*'}
                % increase timerange
                opts.wndlen = opts.wndlen*1.1;                
            case {'divide','/'}
                % decrease timerange
                opts.wndlen = opts.wndlen*0.9;                
            case 'pagedown'
                % shift display page offset down
                opts.pageoffset = opts.pageoffset+1;                
            case 'pageup'
                % shift display page offset up
                opts.pageoffset = opts.pageoffset-1;
            case 'n'
                opts.display_mode = 'new';
            case 'o'
                opts.display_mode = 'old';
            case 'b'
                opts.display_mode = 'both';
            case 'd'
                opts.display_mode = 'diff';
            case 'l'
                opts.show_eventlegend = ~opts.show_eventlegend;
            case 'e'
                opts.show_events = ~opts.show_events;
            case 'h'
                if strcmp(get(hSlider,'Visible'),'on')
                    set(hSlider,'Visible','off')
                else
                    set(hSlider,'Visible','on')
                end
        end        
        on_update();
    end
end

% create a mapping from event types onto colors
function map = gen_colormap(eventstruct,mapname)
if isempty(eventstruct)
    map = struct('keys',[],'values',[]);
else
	map.keys = unique({eventstruct.type});
	if isscalar(map.keys)
		tmp = colormap(mapname);
		map.values = tmp(round(end/2),:);
	elseif ~isempty(map.keys)
		tmp = colormap(mapname);
		map.values = tmp(1+floor((0:length(map.keys)-1)/(length(map.keys)-1)*(length(tmp)-1)),:);
	else
		map.values = [];
	end
end
end
