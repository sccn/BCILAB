function [J,Msp,Nsp]=mtfftpt(data,tapers,nfft,t,f,findx)
% Multi-taper fourier transform for point process given as times
%
% Usage:
% [J,Msp,Nsp]=mtfftpt (data,tapers,nfft,t,f,findx) - all arguments required
% Input: 
%       data        (struct array of times with dimension channels/trials; 
%                   also takes in 1d array of spike times as a column vector) 
%       tapers      (precalculated tapers from dpss) 
%       nfft        (length of padded data) 
%       t           (time points at which tapers are calculated)
%       f           (frequencies of evaluation)
%       findx       (index corresponding to frequencies f) 
% Output:
%       J (fft in form frequency index x taper index x channels/trials)
%       Msp (number of spikes per sample in each channel)
%       Nsp (number of spikes in each channel)
if nargin < 6; error('Need all input arguments'); end;
if isstruct(data); C=length(data); else C=1; end% number of channels
K=size(tapers,2); % number of tapers
nfreq=length(f); % number of frequencies
if nfreq~=length(findx); error('frequency information (last two arguments) inconsistent'); end;
H=fft(tapers,nfft,1);  % fft of tapers
H=H(findx,:); % restrict fft of tapers to required frequencies
w=2*pi*f; % angular frequencies at which ft is to be evaluated
Nsp=zeros(1,C); Msp=zeros(1,C);
for ch=1:C;
  if isstruct(data);
     fnames=fieldnames(data);
     eval(['dtmp=data(ch).' fnames{1} ';'])
     indx=find(dtmp>=min(t)&dtmp<=max(t));
     if ~isempty(indx); dtmp=dtmp(indx);
     end;
  else
     dtmp=data;
     indx=find(dtmp>=min(t)&dtmp<=max(t));
     if ~isempty(indx); dtmp=dtmp(indx);
     end;
  end;
  Nsp(ch)=length(dtmp);
  Msp(ch)=Nsp(ch)/length(t);
  if Msp(ch)~=0;
      data_proj=interp1(t',tapers,dtmp);
      exponential=exp(-i*w'*(dtmp-t(1))');
      J(:,:,ch)=exponential*data_proj-H*Msp(ch);
  else
      J(1:nfreq,1:K,ch)=0;
  end;
end;
