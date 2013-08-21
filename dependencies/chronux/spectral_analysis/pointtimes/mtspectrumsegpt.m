function [S,f,R,varS,zerosp,C,Serr]=mtspectrumsegpt(data,win,params,segave,fscorr)
% Multi-taper segmented spectrum for a univariate binned point process
%
% Usage:
%
% [S,f,R,varS,zerosp,C,Serr]=mtspectrumsegpt(data,win,params,segave,fscorr)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (structure array of one channel of spike times; 
%            also accepts 1d vector of spike times) -- required
%       win  (duration of the segments) - required. 
%       params: structure with fields tapers, pad, Fs, fpass, err
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
%       segave - (0 for don't average over segments, 1 for average) - optional - default  1
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 
%                1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
% Output:
%       S       (spectrum in form frequency x segments if segave=0; function of frequency if segave=1)
%       f       (frequencies)
%       R       (spike rate)
%       varS    (variance of the spectrum as a function of frequency)
%       zerosp  (0 for segments in which spikes were found, 1 for segments
%       C       (covariance matrix of the log spectrum - frequency x
%       frequency matrix)
%       Serr    (error bars) - only if err(1)>=1

if nargin < 2; error('Need data and segment information'); end;
if nargin < 3; params=[]; end;
if nargin < 4 || isempty(segave); segave=1; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear params trialave
if nargin < 5 || isempty(fscorr); fscorr=0;end;

if nargout > 4 && err(1)==0; error('cannot compute error bars with err(1)=0; change params and run again'); end;

dtmp=change_row_to_column(data);
T=max(dtmp); % total length of data
minT=min(dtmp); 
E=minT:win:T-win; % fictitious event triggers
win=[0 win]; % use window length to define left and right limits of windows around triggers
dtmp=createdatamatpt(dtmp,E,win); % create segmented data set
[mintime,maxtime]=minmaxsptimes(dtmp);
dt=1/Fs; % sampling time
t=mintime-dt:dt:maxtime+dt; % time grid for prolates
N=length(t); % number of points in grid for dpss
nfft=max(2^(nextpow2(N)+pad),N); % number of points in fft of prolates
[f,findx]=getfgrid(Fs,nfft,fpass); % get frequency grid for evaluation
tapers=dpsschk(tapers,N,Fs); % check tapers
[J,Msp,Nsp]=mtfftpt(dtmp,tapers,nfft,t,f,findx);% mt fft for point process times
R=Msp*Fs;
S=squeeze(mean(conj(J).*J,2)); % spectra of non-overlapping segments (averaged over tapers)
if segave==1; SS=squeeze(mean(S,2));R=mean(R);else;SS=S;end;% mean of the spectrum averaged across segments
if nargout > 3
    lS=log(SS); % log spectrum for nonoverlapping segments
%     varS=var(lS,1,2); % variance of log spectrum
    varS=var(lS',1)';% variance of the log spectrum R13
    if nargout > 4
       zerosp=zeros(1,size(data,2));
       zerosp(Nsp==0)=1;
       if nargout > 5
          C=cov(lS'); % covariance matrix of the log spectrum
          if nargout==7; 
             if fscorr==1;
                Serr=specerr(SS,J,err,segave,Nsp);
             else
                Serr=specerr(SS,J,err,segave);
             end;
          end;
       end;
    end;
end;
S=SS;
