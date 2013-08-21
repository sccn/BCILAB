function callback_slider_prob(sliderHandle,~,sRate,epoched,epoch_length,ax1,latency,events2plot,events,axLikelihood,windowInS,ylim1,ylim2)

figh = gcf;
currentFrame = sliderHandle.Value;
data = get(figh,'userdata');
windowlength = data{2};
data(2) = [];
data = data{1};
if gca ~= ax1
    axes(ax1);
end

%children = get(ax,'children');
if ~epoched
    
    
    T = ceil(windowlength);
    if sRate*(currentFrame + T)+1 > size(data,2)
        plot(ax1,data(end,(sRate*currentFrame+1:end)),data(1:end-2,(sRate*currentFrame+1:end)))
        set(ax1,'xlim',[currentFrame ceil(size(data,2)/sRate)]);
        xlim = get(ax1,'xlim');
        set(ax1,'xtick',linspace(xlim(1),xlim(2),11))
        
    else
        plot(ax1,data(end,(sRate*currentFrame:sRate*(currentFrame + T))+1),data(1:end-2,(sRate*currentFrame:sRate*(currentFrame + T))+1))
        set(ax1,'xlim',[currentFrame currentFrame+T]);
        set(ax1,'xtick',linspace(currentFrame,currentFrame+T,11));
    end
    
    xlim = get(ax1,'xlim');
    set(ax1,'ylim',[-0.1 1.1]);
    
    ylabel(ax1,'Probability of Model Being Active')
    xlabel(ax1,fastif(epoched,'Trials','Time(sec)'))
    hold(ax1,'on');
    %draw events
    x_coord = [];
    if ~isempty(events2plot)
        count = 1;
        
        for i = 1:length(latency)
            if latency(i)>sRate*xlim(2)
                break;
            else
                if latency(i)>=currentFrame*sRate && ismember({events(i).type},events2plot)
                    x_coord(count) = latency(i)/sRate;
                    plot(ax1,[x_coord(count) x_coord(count)],[-0.1 1.1],'linestyle','- -')
                    ht = text(x_coord(count),1.12,num2str(events(i).type));
                    set(ht,'rotation',90);
                    count = count + 1;
                end
            end
        end
    end
    hold(ax1,'off')
    
    if sRate*(currentFrame + T)+1 > size(data,2)
        plot(axLikelihood,data(end,(sRate*currentFrame+1:end)),data(end-1,(sRate*currentFrame+1:end)))
        set(axLikelihood,'xlim',[currentFrame ceil(size(data,2)/sRate)]);
        
        set(axLikelihood,'xtick',linspace(xlim(1),xlim(2),11))
    else
        plot(axLikelihood,data(end,(sRate*currentFrame:sRate*(currentFrame + T))+1),data(end-1,(sRate*currentFrame:sRate*(currentFrame + T))+1));
        set(axLikelihood,'xlim',[currentFrame currentFrame+T]);
        set(axLikelihood,'xtick',linspace(currentFrame,currentFrame+T,11));
    end
    
    set(axLikelihood,'ylim',[ylim1 ylim2]);
    
    
    xlabel(axLikelihood,fastif(epoched,'Trials','Time(sec)'))
    hold(axLikelihood,'on');
    
    %draw events
    for i = 1:length(x_coord)
        plot(axLikelihood,[x_coord(i) x_coord(i)],[ylim1 ylim2],'linestyle','- -')
    end
    
    hold(axLikelihood,'off')
    if T ~= windowlength
        set(findobj('parent',gcf,'tag','WindowSizeTag'),'string',num2str(T))
    end
    
    sliderHandle.UnitIncrement = ceil(T/10);
    sliderHandle.Maximum = size(data,2)/sRate-T + ceil(T/10);
    sliderHandle.VisibleAmount = sliderHandle.UnitIncrement;
  
    
else
    
    T = ceil(windowlength);
    
    if epoch_length*(currentFrame + T) > size(data,2)
        currentFrame = size(data,2)/epoch_length - T;
        sliderHandle.Value = currentFrame;
    end
    
    plot(ax1,data(end,(currentFrame*epoch_length+1:epoch_length*(currentFrame + T))),data(1:end-2,(currentFrame*epoch_length+1:epoch_length*(currentFrame + T))))
    hold(ax1,'on');
    set(ax1,'ylim',[-0.1 1.1]);
    set(ax1,'xlim',[currentFrame+0.5 currentFrame+T+0.5]);
    lim = get(ax1,'xlim');
    set(ax1,'xtick',floor(lim(1):lim(2)))
    ylabel(ax1,'Probability of Model Being Active')
    xlabel(ax1,fastif(epoched,'Trials','Time(sec)'))
    
    %draw epoch seperating lines
    for i = currentFrame+1:currentFrame+T
        plot(ax1,[i+0.5 i+0.5],[-0.1 1.1],'color',[1 0.5 0.25],'linestyle','- -');
    end
    %-------------------------------
    windowInMs = 1000*windowInS;
    x_coord = [];
    if ~isempty(events2plot)
    % draw event lines
    
    C = 1/(windowInMs(2)-windowInMs(1));
    count = 1;
    
    for i = currentFrame+1:currentFrame+T
        lengthevents = length(events(i).event);
        if lengthevents == 1
            try 
                temp = events(i).eventtype;
            catch
                temp = events(i).eventtype{1};
            end
        else
            temp = '';
        end
       
        for j = 1:lengthevents
            
            if lengthevents == 1
                str = temp{1};
            else
                str = events(i).eventtype{j};
            end
            
            if ismember(num2str(str),events2plot)
            x_coord(count) = C*([events(i).eventlatency{j}]-windowInMs(1)) + (i - 0.5);
            plot(ax1,[x_coord(count) x_coord(count)],[-0.1 1.1],'linestyle','- -')
            ht = text(x_coord(count),1.12,num2str(events(i).eventtype{j}));
            set(ht,'rotation',90);
            count = count + 1;
            end
        end
    end
    
    end
    hold(ax1,'off');
    
    plot(axLikelihood,data(end,(currentFrame*epoch_length+1:epoch_length*(currentFrame + T))),data(end-1,(currentFrame*epoch_length+1:epoch_length*(currentFrame + T))))
    set(axLikelihood,'ylim',[ylim1 ylim2]);
    set(axLikelihood,'xlim',[currentFrame+0.5 currentFrame+T+0.5]);
    set(axLikelihood,'xtick',floor(lim(1):lim(2)))
    xlabel(axLikelihood,fastif(epoched,'Trials','Time(sec)'))
    hold(axLikelihood,'on');
    
    for i = currentFrame+1:currentFrame+T
        plot(axLikelihood,[i+0.5 i+0.5],[ylim1 ylim2],'color',[1 0.5 0.25],'linestyle','- -');
    end
    
    % draw event lines
    
    for i = 1:length(x_coord)
        plot(axLikelihood,[x_coord(i) x_coord(i)],[ylim1 ylim2],'linestyle','- -');
    end
    
    hold(axLikelihood,'off');
    if T ~= windowlength
        set(findobj('parent',gcf,'tag','WindowSizeTag'),'string',num2str(T))
    end
    
    sliderHandle.UnitIncrement = ceil(T/5);
    sliderHandle.Maximum = size(data,2)/epoch_length-T + ceil(T/5);
    sliderHandle.VisibleAmount = sliderHandle.UnitIncrement;
    
end
ylabel(axLikelihood,'Log-likelihood of data under most probable model')






end