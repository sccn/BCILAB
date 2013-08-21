function scrollplot(dx,x,y, epoched,bigploth)
% Created by Steven Lord, slord@mathworks.com
% Uploaded to MATLAB Central
% http://www.mathworks.com/matlabcentral
% 7 May 2002
%
% Permission is granted to adapt this code for your own use.
% However, if it is reposted this message must be intact.

%% dx is the width of the axis 'window'
plot(x,y);
ylabel('Probability of Being Active')
ylim1 = -0.1;
ylim2 = 1.1;
xmax=max(x);
xmin = min(x);
%cb = ['set(gca,''xlim'',get(h,''value'')+[' num2str(xmin) ' ' num2str(xmin+dx) '])'];

if ~epoched
    
    xlabel('Time (sec)')
    
    a=gca;
    
    % Set appropriate axis limits and settings
    set(gcf,'doublebuffer','on');
    %% This avoids flickering when updating the axis
    set(a,'xlim',[0 dx]);
    %set(a,'ylim',[min(y) max(y)]);
    
    % Generate constants for use in uicontrol initialization
    pos=get(a,'position');
    set(a,'position',[pos(1) pos(2)+0.1 pos(3) pos(4)-0.1]);
    %% This will create a slider which is just underneath the axis
    %% but still leaves room for the axis labels above the slider
    
    Newpos=[pos(1) pos(2)-0.05 pos(3) 0.05];
    %S= ['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
    set(gca,'ylim',[ylim1 ylim2]);
    %% Setting up callback string to modify XLim of axis (gca)
    %% based on the position of the slider (gcbo)
    axes(bigploth);
    th = text(0, 1.14,'.','FontSize',20);
    axes(a);
    % Creating Uicontrol
    h=uicontrol('style','slider',...
        'units','normalized','position',Newpos,...
        'min',0,'max',xmax-dx);
    % hJScrollBar = findjobj(h);
    % hJScrollBar.AdjustmentValueChangedCallback = @S;
    set(a,'xlim',[xmin xmin+dx]);
    set(h,'value',1)
    
    jScrollBar = findjobj(h);
    step = jScrollBar.Value + 100000;
    jScrollBar.AdjustmentValueChangedCallback = {@myfunc,dx,xmin,step,a,th};
    
    jScrollBar.Value = -100000;
    
    %lh = addlistener(h,'Action',@myfunc);
    
else
    xlabel('Trials')
    a=gca;
    
    % Set appropriate axis limits and settings
    set(gcf,'doublebuffer','on');
    %% This avoids flickering when updating the axis
    set(a,'xlim',[min(x) min(x)+dx]);
    %set(a,'ylim',[min(y) max(y)]);
    
    % Generate constants for use in uicontrol initialization
    pos=get(a,'position');
    set(a,'position',[pos(1) pos(2)+0.1 pos(3) pos(4)-0.1]);
    %% This will create a slider which is just underneath the axis
    %% but still leaves room for the axis labels above the slider
    
    Newpos=[pos(1) pos(2)-0.05 pos(3) 0.05];
    %S= ['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
    
    set(gca,'ylim',[ylim1 ylim2]);
    %% Setting up callback string to modify XLim of axis (gca)
    %% based on the position of the slider (gcbo)
    axes(bigploth);
    th = text(0, 1.14,'.','FontSize',20);
    axes(a);
    
    for i = 1:floor(xmax)-1
        hl = line([i+min(x) i+min(x)],[ylim1 ylim2]);
        set(hl,'LineStyle','- -');
    end
    % Creating Uicontrol
    h=uicontrol('style','slider',...
        'units','normalized','position',Newpos,...
        'min',0,'max',xmax-dx);
    % hJScrollBar = findjobj(h);
    % hJScrollBar.AdjustmentValueChangedCallback = @S;
    set(gca,'xlim',[xmin xmin+dx]);
    set(gca,'XTick',1:floor(xmax))
    set(h,'value',1)
    jScrollBar = findjobj(h);
    step = jScrollBar.Value + 100000;
    jScrollBar.AdjustmentValueChangedCallback = {@myfunc,dx,xmin,step,a,th};
    jScrollBar.Value = -100000;
    
    %lh = addlistener(h,'Action',@myfunc);
end
end
function myfunc(h,e,dx,xmin,step,a,th, ylim1,ylim2)

set(a,'xlim',(h.Value + 100000)/step+[xmin xmin+dx]);
% if nargin>6
%   lines = get(a,'Children');
%     if length(lines)>3
%        delete(lines(1:dx-1))
%     end
%     for i= 1:dx-1
%     xmin = xmin + (h.Value + 100000)/step;
%     hl = line([i+xmin i+xmin],[ylim1 ylim2]);
%     set(hl,'LineStyle','- -');
%     end
% end
set(th,'Position',[(h.Value + 100000)/step + xmin 1.14])
end


