function [v2plot llt2plot]= modprobplot(EEG, models2plot,smoothlength,events2plot)

if ~isfield(EEG.etc,'amica') || isempty(EEG.etc.amica)
    error('No AMICA solution found. You should first load AMICA components');
end
% if ~isfield(EEG.etc.amica,'v_smooth')
%     error('No AMICA solution found. You should first load AMICA components through Tools>AMICA>Load AMICA components');
% end
if max(models2plot)>EEG.etc.amica.num_models || min(models2plot) < 1
    error('Not a valid range for models')
end

if smoothlength ~= 0
    if isfield(EEG.etc.amica,'smooth_length') && isequal(smoothlength,EEG.etc.amica.smooth_length)
        if isfield(EEG.etc.amica,'v_smooth') && isfield(EEG.etc.amica,'LLt_smooth')
            v2plot = EEG.etc.amica.v_smooth;
            llt2plot = EEG.etc.amica.LLt_smooth;
        else
            disp('Smoothing probabilities...')
            [v2plot llt2plot] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,smoothlength);
        end
    else
        disp('Smoothing probabilities...')
        [v2plot llt2plot] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,smoothlength);
        
    end
else
    smoothlength = [];
    v2plot = EEG.etc.amica.v;
    llt2plot = EEG.etc.amica.LLt;
end

y = v2plot(models2plot,:,:);
figure;
set(gcf,'position',[21 646 560 420]);
if size(v2plot,3) > 1
        T = 5; %number of epochs to be shown per window for scrolling plot.
        size2 = size(y,2);
        size3 = size(y,3);
        t = 1:size2*size3;
        t = (t + size2-ceil(size2/2))/size2;
        plot(t,reshape(y,size(y,1),size(y,3)*size(y,2)))
        h = gca;
        set(h,'ylim',[-0.1 1.1]);
        set(h,'XTick',ceil(linspace(1,max(t),10)))
        xlabel('Trials')
        
            

    
else
        T = 20;% length of the window for scrolling plot.
        t = 1:size(v2plot,2);
        t = (t - 1)/EEG.srate;
        plot(t,y);
        h = gca;
        set(h,'ylim',[-0.1 1.1]);
        %set(h,'XTick',ceil(linspace(1,max(t),10)))
        xlabel('Time (sec)')
        
        
       
    
end
BACKCOLOR = [.93 .96 1];
set(gcf,'color',BACKCOLOR);
ylabel('Probability of Model Being Active')
set(gcf,'name','AMICA model probabilities')
str = {};
for i = 1:length(models2plot)
    newstr = ['Model ' num2str(models2plot(i))];
    str{i} = newstr;
end
legend(str,'Location','NorthEastOutside');
L = max(llt2plot(models2plot,:,:),[],1);
events2plot = cellfun(@num2str,fastif(~iscell(events2plot),{events2plot},events2plot),'UniformOutput',0);
scrollplot_fast(y,L,models2plot,T,EEG.srate, fastif(size(y,3) > 1,EEG.epoch,EEG.event), [EEG.xmin EEG.xmax], events2plot,smoothlength);

%legend(str,'Location','NorthEastOutside');








    
    
