function [SS,tau]=mtspectrum_of_spectrumc(data,win,tapers2spec,params)
% Multi-taper segmented, second spectrum (spectrum of the log spectrum) for a continuous process
% This routine computes the second spectrum by explicitly evaluating the
% Fourier transform (since the spectrum is symmetric in frequency, it uses
% a cosine transform)
%
% Usage:
%
% [SS,tau]=mtspectrum_of_spectrumc(data,win,tapers2spec,params)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (single channel) -- required
%       win  (duration of the segments) - required. 
%       tapers2spec (tapers used for the spectrum of spectrum computation) -
%       required in the form [use TW K] - Note that spectrum of the
%       spectrum involves computing two Fourier transforms. While the first
%       transform (of the original data) is always computed using the
%       multi-taper method, the current routine allows the user to specify 
%       whether or not to use this method for the second transform. use=1
%       means use tapers, use=anything other than 1 means do not use the
%       multitaper method. If use=1, then tapers2spec controls the
%       smoothing for the second Fourier transform. Otherwise, a direct
%       Fourier transform is computed.
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
%                                   Default all frequencies between 0 and
%                                   Fs/2
% Output:
%       SS       (second spectrum in form frequency x segments x trials x channels 
%                if segave=0; in the form frequency x trials x channels if segave=1)
%       tau      (frequencies)
if nargin < 3; error('Need data,segment duration and taper information'); end;
if nargin < 4 ; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params); 
[N,Ntr,NC]=size(data);
if Ntr==1; error('cannot compute second spectrum with just one trial'); end;
dt=1/Fs; % sampling interval
T=N*dt; % length of data in seconds
E=0:win:T-win; % fictitious event triggers
datatmp=createdatamatc(data(:,1,1),E,Fs,[0 win]); % segmented data
Ninseg=size(datatmp,1); % number of samples in segments
nfft=max(2^(nextpow2(Ninseg)+pad),Ninseg);
[f,findx]=getfgrid(Fs,nfft,fpass); 
NF=length(findx);
S=zeros(NF,Ntr,NC);
for nc=1:NC;
    for ntr=1:Ntr;
        datatmp=change_row_to_column(data(:,ntr,nc));
        s=mtspectrumsegc(datatmp,win,params,1);
        S(:,ntr,nc)=s;
    end
end;
Sm=mean(S,2);
if use==1;
   params.tapers=tapers2spec;
   params.Fs=1/(f(2)-f(1));
   params.fpass=[0 params.Fs/2];
else;
   tau=[0:NF-1]/max(f);
   cosinefunc=cos(2*pi*f'*tau);
end;

for nc=1:NC;
    for ntr=1:Ntr;
        s=S(:,ntr,nc)./Sm(:,nc);
        s=log(s);
        if use==1;
            sflip=flipdim(s,1);
            s=[sflip(1:NF-1);s];
            [ss,tau]=mtspectrumc(s,params);
            SS(:,ntr,nc)=ss;
        else;
            s=repmat(s,[1 NF]).*cosinefunc;
    %         subplot(221); plot(s(:,1));
    %         subplot(222); plot(s(:,10));
    %         subplot(223); plot(s(:,100));
    %         subplot(224); plot(s(:,120));
    %         pause
            s=trapz(f,s,1)';
            ss=s.*conj(s);
%         plot(tau,s)
%         pause
        end
        SS(:,ntr,nc)=ss;
    end
end;
SS=mean(SS,2);
