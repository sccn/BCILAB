function [P,Fstat,f0]=mtpowerandfstatc(data,params,f0)
% Multi-taper computation of the power and the fstatistic for a particular frequency - continuous process
%
% Usage:
%
% [P,Fstat,f0]=mtpowerandfstatc(data,params,f0)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (in form samples x channels/trials or a single vector) -- required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
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
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%       f0  (frequency of calculation)
% Output:
%       P       (integrated power within the frequency range of interest (trapezoidal integration))
%       Fstat   (F-statistic)
%       f0      (frequency)

if nargin < 1; error('Need data'); end;
if nargin < 2; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear fpass err trialave params
data=change_row_to_column(data);
[N,C]=size(data);
tapers=dpsschk(tapers,N,Fs); % calculate the tapers
[N,K]=size(tapers);
nfft=max(2^(nextpow2(N)+pad),N);% number of points in fft
%[f0,findx]=getfgrid(Fs,nfft,f0);% frequency grid to be returned

tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
data=data(:,:,ones(1,K)); % add taper indices to data
data=permute(data,[1 3 2]); % reshape data to get dimensions to match those of tapers
data_proj=data.*tapers; % product of data with tapers in the form time x tapers x channels
t=(0:N-1)'/Fs;
fourier=exp(-i*2*pi*f0*t);
fourier=fourier(:,ones(1,K),ones(1,C));
J=squeeze(sum(fourier.*data_proj))/Fs; 

Kodd=1:2:K;
Keven=2:2:K;
tapers=tapers(:,:,ones(1,C)); % add channel indices to the tapers - t x K x C
H0 = squeeze(sum(tapers(:,Kodd,:),1)); % calculate sum of tapers for even prolates - K x C 

if C==1; H0=H0'; J=J'; end;
P=squeeze(mean(J.*conj(J),1));
Jp=J(Kodd,:); % drop the even ffts
H0sq=sum(H0.*H0,1);% sum of squares of H0^2 across taper indices - dimensions C
JpH0=sum(Jp.*H0,1);% sum of the product of Jp and H0 across taper indices - f x C\
A=squeeze(JpH0./H0sq); % amplitudes for all frequencies and channels
Kp=size(Jp,1); % number of even prolates
Ap=A(ones(1,Kp),:); % add the taper index to C
Jhat=Ap.*H0; % fitted value for the fft

num=(K-1).*(abs(A).^2).*squeeze(H0sq);%numerator for F-statistic
den=squeeze(sum(abs(Jp-Jhat).^2,1)+sum(abs(J(Keven,:)).^2,1));% denominator for F-statistic
Fstat=num./den; % F-statisitic

