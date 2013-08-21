function [S,t,f,R,Serr]=mtspecgrampb(data,movingwin,params,fscorr)
% Multi-taper time-frequency spectrum - binned point process
%
% Usage:
%
% [S,t,f,R,Serr]=mtspecgrampb(data,movingwin,params,fscorr)
% Input: 
%       data        (in form samples x channels/trials or single vector) -- required
%       movingwin         (in the form [window,winstep] i.e length of moving
%                                                 window and step size.
%                                                 
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms: 
%                   (1) A numeric vector [TW K] where TW is the
%                       time-bandwidth product and K is the number of
%                       tapers to be used (less than or equal to
%                       2TW-1). 
%                   (2) A numeric vector [W T p] where W is the
%                       bandwidth, T is the duration of the data and p 
%                       is an integer such that 2TW-p tapers are used. In
%                       this form there is no default i.e. to specify
%                       the bandwidth, you have to specify T and p as
%                       well. Note that the units of W and T have to be
%                       consistent: if W is in Hz, T must be in seconds
%                       and vice versa. Note that these units must also
%                       be consistent with the units of params.Fs: W can
%                       be in Hz if and only if params.Fs is in Hz.
%                       The default is to use form 1 with TW=3 and K=5
%                   Note that T has to be equal to movingwin(1).
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials/channnels when 1, don't average when 0) - optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 
%                1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
% Output:
%       S       (spectrum in form time x frequency x channels/trials for trialave=0; 
%               or as a function of frequency if trialave=1)
%       t       (times)
%       f       (frequencies)
%       R       (rate)
%       Serr    (error bars) - only for err(1)>=1

if nargin < 2; error('Need data and window parameters'); end;
if nargin < 3; params=[]; end;

if length(params.tapers)==3 & movingwin(1)~=params.tapers(2);
    error('Duration of data in params.tapers is inconsistent with movingwin(1), modify params.tapers(2) to proceed')
end

[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);

if nargin < 4 || isempty(fscorr); fscorr=0; end;
if nargout > 4 && err(1)==0; 
%    error('Cannot compute errors with err(1)=0'); 
     error('When Serr is desired, err(1) has to be non-zero.');
end;
data=change_row_to_column(data);
[N,Ch]=size(data);
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=max(2^(nextpow2(Nwin)+pad),Nwin);
f=getfgrid(Fs,nfft,fpass); Nf=length(f);
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers
winstart=1:Nstep:N-Nwin+1;
nw=length(winstart);
if trialave
    S = zeros(nw,Nf);
    R = zeros(nw,Nwin);
    if nargout==4; Serr=zeros(2,nw,Nf); end;
else
    S = zeros(nw,Nf,Ch);
    R = zeros(nw,Nwin,Ch);
    if nargout==4; Serr=zeros(2,nw,Nf,Ch); end;
end

for n=1:nw;
   indx=winstart(n):winstart(n)+Nwin-1;
   datawin=data(indx,:);
   if nargout==5;
     [s,f,r,serr]=mtspectrumpb(datawin,params,fscorr);
     Serr(1,n,:,:)=squeeze(serr(1,:,:));
     Serr(2,n,:,:)=squeeze(serr(2,:,:));
   else
     [s,f,r]=mtspectrumpb(datawin,params,fscorr);
   end;
   S(n,:,:)=s;
   R(n,:)=r';
end;
winmid=winstart+round(Nwin/2);
t=winmid/Fs;
S=squeeze(S); R=squeeze(R); if nargout==5; Serr=squeeze(Serr);end
