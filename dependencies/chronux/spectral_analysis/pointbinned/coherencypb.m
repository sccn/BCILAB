function [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencypb(data1,data2,params,fscorr)
% Multi-taper coherency,cross-spectrum and individual spectra - binned point process
%
% Usage:
%
% [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencypb(data1,data2,params,fscorr)
% Input: 
%       data1 (in form samples x trials) -- required
%       data2 (in form samples x trials) -- required
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
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 
%                 1 (use finite size corrections) - optional  
%                 (available only for spikes). Defaults 0.
% Output:
%       C (magnitude of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       phi (phase of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S12 (cross spectrum -  frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S1 (spectrum 1 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S2 (spectrum 2 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       f (frequencies)     
%       zerosp (1 for trials in either channel where spikes were absent, zero otherwise)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phistd - jackknife/theoretical standard deviation for phi.  Note that 
%                phi + 2 phistd and phi -2 phistd will give 95% confidence bands 
%                for phi - only for err(1)>=1
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)
if nargin < 2; error('Need data1 and data2'); end;
if nargin < 2; error('Need data1 and data2'); end;
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear params
if nargin < 4 || isempty(fscorr); fscorr=0; end;

if nargout > 9 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs or outputs and run again');
end;
if nargout > 7 && err(1)==0;
    error('When errors are desired, err(1) has to be non-zero.');
end;

[N,Ch]=check_consistency(data1,data2);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
[J1,Msp,Nsp1]=mtfftpb(data1,tapers,nfft);
[J2,Msp,Nsp2]=mtfftpb(data2,tapers,nfft);
zerosp=zeros(1,Ch); % initialize the zerosp variable
zerosp(Nsp1==0 | Nsp2==0)=1; % set the zerosp variable
J1=J1(findx,:,:);
J2=J2(findx,:,:);
S12=squeeze(mean(conj(J1).*J2,2));
S1=squeeze(mean(conj(J1).*J1,2));
S2=squeeze(mean(conj(J2).*J2,2));
if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;
C12=S12./sqrt(S1.*S2);
C=abs(C12);
phi=angle(C12);
if nargout==10; 
    if fscorr==1; 
       [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave,Nsp1,Nsp2);
    else
       [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave);
    end;
elseif nargout==9;
    if fscorr==1; 
        [confC,phistd]=coherr(C,J1,J2,err,trialave,Nsp1,Nsp2);
    else
        [confC,phistd]=coherr(C,J1,J2,err,trialave);
    end;
end;
clear Msp
