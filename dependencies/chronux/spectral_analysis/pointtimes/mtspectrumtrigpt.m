function [S,f,R,Serr]=mtspectrumtrigpt(data,E,win,params,fscorr)
% Multi-taper time-frequency spectrum - point process times
%
% Usage:
%
% [S,f,R,Serr]=mtspectrumtrigpt(data,E,win,params,fscorr)
% Input: 
%       data        (structure array of one channel of spike times; 
%                   also accepts 1d column vector of spike times) -- required
%       E           (event times) - required
%       win         (in the form [winl winr] i.e window around each event)--
%                                                 required
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
%           trialave (average over events when 1, don't average when 0) -
%           default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 
%                1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
%
% Output:
%       S       (triggered spectrum in form frequency x events if trialave=0; 
%               function of frequency if trialave=1)
%       f       (frequencies)
%       R       (spike rate)
%       Serr    (error bars) - only for err(1)>=1

if nargin < 3; error('Need data, events and window parameters'); end;
if nargin < 2; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear tapers pad Fs fpass trialave
if nargin < 5 || isempty(fscorr); fscorr=0; end;
if nargout > 3 && err(1)==0; error('Cannot compute errors if err(1)=0'); end;
data=createdatamatpt(data,E,win);
if nargout==4; 
    [S,f,R,Serr]=mtspectrumpt(data,params,fscorr);
else
    [S,f,R]=mtspectrumpt(data,params,fscorr);
end;
