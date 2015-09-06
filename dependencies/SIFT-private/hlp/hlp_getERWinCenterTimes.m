function wintimes = hlp_getERWinCenterTimes(EEG,params)
% return the timepoints (seconds) corresponding to the centers of windows
% w.r.t stimulus event

startp  = 1;
endp    = EEG.CAT.pnts;
winlen  = floor(params.winlen*EEG(1).srate);
winstep = params.winstep*EEG(1).srate;

[dummy eventT] = getindex(EEG(1).times,0);
winpoints      = floor((startp:winstep:endp-winlen)*EEG.srate)+1;
% winpoints =startp:winstep:endp-winlen+1;
wintimes       = (winpoints-eventT+winlen/2)/EEG.srate;  % centers of windows