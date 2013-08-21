function [dS,f]=mtdspectrumpb(data,phi,params)
% Multi-taper spectral derivative - binned point process
%
% Usage:
%
% [dS,f]=mtdspectrumpb(data,phi,params)
% Input: 
%   Note that all times can be in arbitrary units. But the units have to be
%   consistent. So, if E is in secs, win, t have to be in secs, and Fs has to
%   be Hz. If E is in samples, so are win and t, and Fs=1. In case of spike
%   times, the units have to be consistent with the units of data as well.
%       data (in form samples x channels/trials or single vector) -- required
%       tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. 
%                                                 If not specified, use [NW K]=[3 5]
%       phi         (angle for evaluation of derivative) -- required.
%                       e.g. phi=[0,pi/2] giving the time and frequency
%                       derivatives
%       params: structure with fields tapers, pad, Fs, fpass, trialave
%       -optional
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
%                                   Default all frequencies between 0 and
%                                   Fs/2
%           trialave (average over trials when 1, don't average when 0) -
%           optional. Default 0
% Output:
%       dS       (derivative of the spectrum in form phi x frequency x channels/trials if trialave=0; 
%                in the form phi x frequency if trialave=1)
%       f        (frequencies)

if nargin < 2; error('Need data and angle'); end;
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear err params
data=change_row_to_column(data);
N=size(data,1);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass);
tapers=dpsschk(tapers,N,Fs); % check tapers
K=size(tapers,2);
J=mtfftpb(data,tapers,nfft);  
J=J(findx,:,:);
A=sqrt(1:K-1);
A=repmat(A,[size(J,1) 1]);
A=repmat(A,[1 1 size(J,3)]);
% S=squeeze(mean(J(:,1:K-1,:).*conj(J(:,2:K,:)),2));
S=squeeze(mean(J(:,1:K-1,:).*A.*conj(J(:,2:K,:)),2));
if trialave; S=squeeze(mean(S,2)); end;
nphi=length(phi);
for p=1:nphi;
    dS(p,:,:)=real(exp(i*phi(p))*S);
end;
dS=squeeze(dS);
dS=change_row_to_column(dS);
