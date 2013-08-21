function [datac,datafit,Amps,freqs]=rmlinesmovingwinc(data,movingwin,tau,params,p,plt,f0)
% fits significant sine waves to data (continuous data) using overlapping windows.
%
% Usage: [datac,datafit]=rmlinesmovingwinc(data,movingwin,tau,params,p,plt)
%
%  Inputs:  
% Note that units of Fs, fpass have to be consistent.
%       data        (data in [N,C] i.e. time x channels/trials or as a single vector) - required.
%       movingwin         (in the form [window winstep] i.e length of moving
%                                                 window and step size)
%                                                 Note that units here have
%                                                 to be consistent with
%                                                 units of Fs - required
%       tau      parameter controlling degree of smoothing for the amplitudes - we use the
%       function 1-1/(1+exp(-tau*(x-Noverlap/2)/Noverlap) in the region of overlap to smooth
%       the sinewaves across the overlap region. Noverlap is the number of points 
%       in the overlap region. Increasing tau leads to greater overlap smoothing, 
%       typically specifying tau~10 or higher is reasonable. tau=1 gives an almost
%       linear smoothing function. tau=100 gives a very steep sigmoidal. The default is tau=10.
%       params      structure containing parameters - params has the
%       following fields: tapers, Fs, fpass, pad
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
%                    Note that T has to be equal to movingwin(1).
%
%	        Fs 	        (sampling frequency) -- optional. Defaults to 1.
%               fpass       (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%	    p		    (P-value to calculate error bars for) - optional.
%	    Defaults to 0.05/Nwin where Nwin is length of window which
%	    corresponds to a false detect probability of approximately 0.05.
%       plt         (y/n for plot and no plot respectively) - default no
%                   plot.
%       f0          frequencies at which you want to remove the
%                   lines - if unspecified the program uses the f statistic
%                   to determine appropriate lines.
%
%  Outputs: 
%       datafit        (fitted sine waves)
%       datac          (cleaned up data)
if nargin < 2; error('Need data and window parameters'); end;
if nargin < 4 || isempty(params); params=[]; end; 

if length(params.tapers)==3 & movingwin(1)~=params.tapers(2);
    error('Duration of data in params.tapers is inconsistent with movingwin(1), modify params.tapers(2) to proceed')
end

[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params); % set defaults for params
clear err trialave
if nargin < 6; plt='n'; end;
%
% Window,overlap and frequency information
%
data=change_row_to_column(data);
[N,C]=size(data);
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
Noverlap=Nwin-Nstep; % number of points in overlap
%
% Sigmoidal smoothing function
%
if nargin < 3 || isempty(tau); tau=10; end; % smoothing parameter for sigmoidal overlap function
x=(1:Noverlap)';
smooth=1./(1+exp(-tau.*(x-Noverlap/2)/Noverlap)); % sigmoidal function
smooth=repmat(smooth,[1 C]);
%
% Start the loop
%
if nargin < 5 || isempty(p); p=0.05/Nwin; end % default for p value
if nargin < 7 || isempty(f0); f0=[]; end; % empty set default for f0 - uses F statistics to determine the frequencies
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers
winstart=1:Nstep:N-Nwin+1;
nw=length(winstart); 
datafit=zeros(winstart(nw)+Nwin-1,C);
Amps=cell(1,nw);
freqs=cell(1,nw);
for n=1:nw;
   indx=winstart(n):winstart(n)+Nwin-1;
   datawin=data(indx,:);
   [datafitwin,as,fs]=fitlinesc(datawin,params,p,'n',f0);
   Amps{n}=as;
   freqs{n}=fs;
   datafitwin0=datafitwin;
   if n>1; datafitwin(1:Noverlap,:)=smooth.*datafitwin(1:Noverlap,:)+(1-smooth).*datafitwin0(Nwin-Noverlap+1:Nwin,:);end;
   datafit(indx,:)=datafitwin;
end;
datac=data(1:size(datafit,1),:)-datafit;     
if strcmp(plt,'y');
    [S,f]=mtspectrumsegc(data,movingwin(1),params);
    [Sc,fc]=mtspectrumsegc(datac,movingwin(1),params);
    plot(f,10*log10(S),fc,10*log10(Sc));
end;