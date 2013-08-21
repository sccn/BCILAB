function [dS,t,f]=mtdspecgramc(data,movingwin,phi,params)
% Multi-taper derivative of the time-frequency spectrum - continuous process
%
% Usage:
%
% [dS,t,f]=mtdspecgramc(data,movingwin,phi,params)
% Input: 
%   Note that all times can be in arbitrary units. But the units have to be
%   consistent. So, if E is in secs, win, t have to be in secs, and Fs has to
%   be Hz. If E is in samples, so are win and t, and Fs=1. In case of spike
%   times, the units have to be consistent with the units of data as well.
%
%       data        (in form samples x channels/trials or a single vector) -- required
%       movingwin         (in the form [window winstep] i.e length of moving
%                                                 window and step size.
%                                                 Note that units here have
%                                                 to be consistent with
%                                                 units of Fs - required
%       phi         (angle for evaluation of derivative) -- required
%                       e.g. phi=[0,pi/2] giving the time and frequency
%                       derivatives
%       params: structure with fields tapers, pad, Fs, fpass, trialave
%       -optional
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms: 
%                    (1) A numeric vector [TW K] where TW is the
%                        time-bandwidth product and K is the number of
%                        tapers to be used (less than or equal to
%                        2TW-1). 
%                    (2) A numeric vector [W T p] where W is the
%                        bandwidth, T is the duration of the data and p 
%                        is an integer such that 2TW-p tapers are used. In
%                        this form there is no default i.e. to specify
%                        the bandwidth, you have to specify T and p as
%                        well. Note that the units of W and T have to be
%                        consistent: if W is in Hz, T must be in seconds
%                        and vice versa. Note that these units must also
%                        be consistent with the units of params.Fs: W can
%                        be in Hz if and only if params.Fs is in Hz.
%                        The default is to use form 1 with TW=3 and K=5
%                     Note that T has to be equal to movingwin(1).
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
%           trialave - (average over trials/channels when 1, don't average when 0) - optional. Default 0
% Output:
%       dS      (spectral derivative in form phi x time x frequency x channels/trials if trialave=0; 
%               in form phi x time x frequency if trialave=1)
%       t       (times)
%       f       (frequencies)

if nargin < 3; error('Need data, window parameters and angle'); end;
if nargin < 4; params=[]; end;

if length(params.tapers)==3 & movingwin(1)~=params.tapers(2);
    error('Duration of data in params.tapers is inconsistent with movingwin(1), modify params.tapers(2) to proceed')
end

[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear err
data=change_row_to_column(data);
[N,C]=size(data);
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=max(2^(nextpow2(Nwin)+pad),Nwin);
f=getfgrid(Fs,nfft,fpass); Nf=length(f);
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers
params.tapers=tapers;
winstart=1:Nstep:N-Nwin+1;
nw=length(winstart);
if trialave==0; dS=zeros(length(phi),nw,Nf,C); else dS=zeros(length(phi),nw,Nf); end; 
for n=1:nw;
   indx=winstart(n):winstart(n)+Nwin-1;
   datawin=data(indx,:);
   [ds,f]=mtdspectrumc(datawin,phi,params);
   dS(:,n,:,:)=ds;
end;
dS=squeeze(dS);
sz=size(dS);
% if length(sz)==3;
%    dS=permute(dS,[2 1 3 4]);
% elseif length(phi)>1
%    dS=permute(dS,[2 1 3]);
% end;
winmid=winstart+round(Nwin/2);
t=winmid/Fs;
