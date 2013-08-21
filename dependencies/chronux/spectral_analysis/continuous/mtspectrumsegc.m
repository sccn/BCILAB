function [S,f,varS,C,Serr]=mtspectrumsegc(data,win,params,segave)
% Multi-taper segmented spectrum for a univariate continuous process
%
% Usage:
%
% [S,f,varS,C,Serr]=mtspectrumsegc(data,win,params,segave)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (single channel) -- required
%       win  (duration of the segments) - required. 
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
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
%           trialave - not used
%       segave - optional 0 for don't average over segments, 1 for average - default
%       1
% Output:
%       S       (spectrum in form frequency x segments if segave=0; in the form frequency if segave=1)
%       f       (frequencies)
%       varS    (variance of the log spectrum)
%       C       (covariance matrix of the log spectrum - frequency x
%       frequency matrix)
%       Serr    (error bars) only for err(1)>=1

if nargin < 2; error('Need data and segment information'); end;
data=change_row_to_column(data);
if size(data,2)~=1; error('works for only univariate time series'); end;
if nargin < 3 ; params=[]; end;
if nargin < 4 || isempty(segave); segave=1; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params); clear trialave params
if nargout==4 && err(1)==0; 
%   Errors can't be computed if err(1)=0. Need to change params and run again.
    error('When Serr is desired, err(1) has to be non-zero.');
end;
N=size(data,1); % length of segmented data
dt=1/Fs; % sampling interval
T=N*dt; % length of data in seconds
E=0:win:T-win; % fictitious event triggers
win=[0 win]; % use window length to define left and right limits of windows around triggers
data=createdatamatc(data,E,Fs,win); % segmented data
N=size(data,1); % length of segmented data
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
J=mtfftc(data,tapers,nfft,Fs); % compute tapered fourier transforms
J=J(findx,:,:); % restrict to specified frequencies
S=squeeze(mean(conj(J).*J,2)); % spectra of non-overlapping segments (average over tapers)
if segave==1; SS=squeeze(mean(S,2)); else; SS=S;end; % mean of the spectrum averaged across segments
if nargout > 2
    lS=log(S); % log spectrum for nonoverlapping segments
    varS=var(lS',1)'; % variance of log spectrum
%     varS=var(lS',1)';% variance of the log spectrum R13
    if nargout > 3
       C=cov(lS'); % covariance matrix of the log spectrum
       if nargout==5; 
          Serr=specerr(SS,J,err,segave);
       end;
    end;
end;
S=SS;
