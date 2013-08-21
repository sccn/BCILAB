function [Sc,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMatc(data,win,params)
% Multi-taper cross-spectral matrix - another routine, allows for multiple trials and channels 
% Does not do confidence intervals. Also this routine always averages over trials - continuous process
%
% Usage:
%
% [Sc,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMatc(data,win,params)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (in form samples x channels x trials) 
%       win  (duration of non-overlapping window)
%       params: structure with fields tapers, pad, Fs, fpass
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
%	        pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	 e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	 to 512 points; if PAD = 2, we pad the FFT
%			      	 to 2048 points, etc.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
% Output:
%       Sc (cross spectral matrix frequency x channels x channels)
%       Cmat Coherence matrix frequency x channels x channels
%       Ctot Total coherence: SV(1)^2/sum(SV^2) (frequency)
%       Cvec leading Eigenvector (frequency x channels)
%       Cent A different measure of total coherence: GM/AM of SV^2s
%       f (frequencies)  
d=ndims(data);
if d<2, error('Need multidimensional array'); end
if d==2, [N,C]=size(data); end;
if d==3, [N,C,Ntr]=size(data); end; 
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear err trialave params
nwin=round(win*Fs); nfft=max(2^(nextpow2(nwin)+pad),nwin); 
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,nwin,Fs); % check tapers
Sc=zeros(length(findx),C,C);

Nwins=floor(N/nwin);

if d==3, % If there are multiple trials
for iwin=1:Nwins,
    for i=1:Ntr, 
        data1=squeeze(data(1+(iwin-1)*nwin:iwin*nwin,:,i));
        J1=mtfftc(detrend(data1),tapers,nfft,Fs);
        J1=J1(findx,:,:);
        for k=1:C,
            for l=1:C,
                spec=squeeze(mean(conj(J1(:,:,k)).*J1(:,:,l),2)); 
            Sc(:,k,l)=Sc(:,k,l)+spec;
            end
        end
    end
end
Sc=Sc/(Nwins*Ntr);
end

if d==2, % only one trial
    if Nwins ==1
        % fast version
        J1=mtfftc(data,tapers,nfft,Fs);
        J1=J1(findx,:,:);
        Sc = squeeze(sum(bsxfun(@times,conj(J1),permute(J1,[1 2 4 3])),2))./size(tapers,2);
    else
        for iwin=1:Nwins,
            data1=squeeze(data(1+(iwin-1)*nwin:iwin*nwin,:));
            J1=mtfftc(data1,tapers,nfft,Fs);
            J1=J1(findx,:,:);
            for k=1:C,
                for l=1:C,
                    Sc(:,k,l)=Sc(:,k,l)+squeeze(mean(conj(J1(:,:,k)).*J1(:,:,l),2));
                end
            end
        end
        Sc=Sc/Nwins;
    end
end

if nargout > 1
    Cmat=Sc;
    Sdiag=zeros(length(findx),C);
    for k=1:C,
        Sdiag(:,k)=squeeze(Sc(:,k,k));
    end

    for k=1:C,
        for l=1:C,
            Cmat(:,k,l)=Sc(:,k,l)./sqrt(abs(Sdiag(:,k).*Sdiag(:,l)));
        end
    end
end

if nargout > 2
    Ctot=zeros(length(findx),1); Cent=Ctot;
    Cvec=zeros(length(findx),C);
    for i=1:length(findx),
        [u s]=svd(squeeze(Sc(i,:,:)));s=diag(s);
        Ctot(i)=s(1).^2/sum(s.^2); Cent(i)=exp(mean(log(s.^2)))/mean(s.^2);
        Cvec(i,:)=transpose(u(:,1));
    end
end
