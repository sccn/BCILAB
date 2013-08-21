[fname,pname]=uigetfile('*');
eval(['load ' pname '\' fname]);
params.tapers=[5 9]; params.pad=2; params.Fs=1000; params.fpass=[0 100]; params.trialave=1;
movingwin=[0.5 0.1]; % duration of moving window used to evaluate spectrograms
win=[2 2]; % window around events
fignum=[2 3 6 9 8 7 4 1]; % figure numbers for particular directions
for targ=0:7;
    E=targoff(find(targets==targ));
    [Ssp,tsp,fsp,R]=mtspecgramtrigpt(dsp(1),E,win,movingwin,params);
    figure(1);
    subplot(3,3,fignum(targ+1));
    imagesc(t,f,10*log10(S)'); axis xy; colorbar;
    figure(2);
    subplot(3,3,fignum(targ+1));
    plot(t,R); 
    [Slfp,tlfp,flfp]=mtspecgramtrigc(dlfp(:,1),E,win,movingwin,params);
    figure(3);
    subplot(3,3,fignum(targ+1));
    imagesc(t,f,10*log10(S)'); axis xy; colorbar;
end;
    
    
    
    
    


