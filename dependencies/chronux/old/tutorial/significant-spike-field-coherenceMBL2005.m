[fname,pname]=uigetfile('*');
eval(['load ' pname '\' fname]);
params.tapers=[5 9]; params.pad=2; params.Fs=1000; params.fpass=[0 100]; params.trialave=1;params.err=[1 0.05];
movingwin=[0.5 0.1]; % duration of moving window used to evaluate spectrograms
win=[2 2]; % window around events
fignum=[2 3 6 9 8 7 4 1]; % figure numbers for particular directions
for targ=0:7;
    E=targoff(find(targets==targ));
    datasp=createdatamatpt(dsp(1),E,win);
    datalfp=createdatamatc(dlfp(:,1),E,Fs,win);
    [C,phi,S12,S1,S2,t,f,zerosp,confC,phierr]=cohgramcpt(datalfp,datasp,movingwin,params); % Note cohgramcpt does give you spectra but not the errors.
    subplot(3,3,fignum(targ+1));plotsig(C,confC,t,f); axis xy; colorbar;
end;
