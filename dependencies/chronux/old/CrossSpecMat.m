function [Sc,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMat(data,win,params)
%
%
% Multi-taper cross-spectral matrix - another routine, this one allows for multiple trials and channels 
% but does not do confidence intervals. Also this routine always averages over trials - continuous process
%
% Usage:
%
% [Sc,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMat(data,win,params)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (in form samples x channels x trials) 
%       win  (duration of non-overlapping window)
%       params: structure with fields tapers, pad, Fs, fpass
%       - optional
%           tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. 
%                                                 If not specified, use [NW K]=[3 5]
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
nwin=round(win*Fs); nfft=2^(nextpow2(nwin)+pad); 
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

Ctot=zeros(length(findx),1); Cent=Ctot;
Cvec=zeros(length(findx),C);
for i=1:length(findx),
    [u s]=svd(squeeze(Sc(i,:,:)));s=diag(s);
    Ctot(i)=s(1).^2/sum(s.^2); Cent(i)=exp(mean(log(s.^2)))/mean(s.^2);             
    Cvec(i,:)=transpose(u(:,1));

end
        
