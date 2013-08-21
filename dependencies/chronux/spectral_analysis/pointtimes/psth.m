function [R,t,E] = psth(data,sig,plt,T,err,t)
%  function to plot trial averaged rate smoothed by 
%  a Gaussian kernel - visual check on stationarity 
%  Usage: [R,t,E] = psth(data,sig,plt,T,err,t)
%                                                   
%  Inputs:                                  
% Note that all times have to be consistent. If data
% is in seconds, so must be sig and t. If data is in 
% samples, so must sig and t. The default is seconds.
% data            structural array of spike times   
% sig             std dev of Gaussian (default 50ms)
%                 (minus indicates adaptive with    
%                  approx width equal to mod sig)   
% plt = 'n'|'r' etc      (default 'r')              
% T is the time interval (default all)              
% err - 0 = none                                    
%       1 = Poisson                                 
%       2 = Bootstrap over trials (default)         
% (both are based on 2* std err rather than 95%)    
% t   = times to evaluate psth at                   
%                                                   
% The adaptive estimate works by first estimating   
% the psth with a fixed kernel width (-sig) it      
% then alters the kernel width so that the number   
% of spikes falling under the kernel is the same on 
% average but is time dependent.  Reagions rich     
% in data therefore have their kernel width reduced 
%                                                    
% Outputs:                                 
%                                                   
% R = rate                                          
% t = times                                         
% E = errors (standard error)                       

if nargin <1 ; error('I need data!');end
[data]=padNaN(data); % create a zero padded data matrix from input structural array
sz=size(data);
% transposes data so that the longer dimension is assumed to be time
% Procedure used to bring input data into form compatible with that required
% for Murray's routine
if sz(1)>sz(2); data=data';end; 
if nargin < 2; sig = 0.05;end
if nargin < 3; plt = 'r';end
if nargin < 5; err = 2; end

if isempty(sig); sig = 0.05;end
if isempty(plt); plt = 'r';end
if isempty(err); err = 2; end

adapt = 0;
if sig < 0; adapt = 1; sig = -sig; end

%  to avoid end effects increase time interval by the width
%  of the kernel (otherwise rate is too low near edges)

if nargin < 4; 
  T(1) = min(data(:,1));
  T(2) = max(max(data));
else
  T(1) = T(1)-4*sig;
  T(2) = T(2)+4*sig;
end

% calculate NT and ND and filter out times of interest

NT = length(data(:,1));
if NT < 4 && err == 2
  disp('Using Poisson errorbars as number of trials is too small for bootstrap')
  err = 1;
end    
    
m = 1;
D = zeros(size(data));
ND=zeros(1,NT);
for n=1:NT
  indx = find(~isnan(data(n,:)) & data(n,:)>=T(1) & data(n,:)<=T(2));
  ND(n) = length(indx);
  D(n,1:ND(n)) = data(n,indx);
  m = m + ND(n); 
end
N_tot = m;
N_max = max(ND);
D = D(:,1:N_max);

% if the kernel density is such that there are on average 
% one or less spikes under the kernel then the units are probably wrong

L = N_tot/(NT*(T(2)-T(1)));
if 2*L*NT*sig < 1 || L < 0.1 
  disp('Spikes very low density: are the units right? is the kernel width sensible?')
  disp(['Total events: ' num2str(fix(100*N_tot)/100) ' sig: ' ...
        num2str(fix(1000*sig)) 'ms T: ' num2str(fix(100*T)/100) ' events/sig: ' ...
        num2str(fix(100*N_tot*sig/(T(2)-T(1)))/100)])
end

%    Smear each spike out  
%    std dev is sqrt(rate*(integral over kernal^2)/trials)     
%    for Gaussian integral over Kernal^2 is 1/(2*sig*srqt(pi))

if nargin < 6
  N_pts =  fix(5*(T(2)-T(1))/sig);
  t = linspace(T(1),T(2),N_pts);
else
  N_pts = length(t);
end
  
RR = zeros(NT,N_pts);
f = 1/(2*sig^2);
for n=1:NT
  for m=1:ND(n)
    RR(n,:) = RR(n,:) + exp(-f*(t-D(n,m)).^2);
  end
end
RR = RR*(1/sqrt(2*pi*sig^2));
if NT > 1; R = mean(RR); else R = RR;end

if err == 1
  E = sqrt(R/(2*NT*sig*sqrt(pi)));
elseif err == 2
  Nboot = 10;
  mE = 0;
  sE = 0;
  for b=1:Nboot
    indx = floor(NT*rand(1,NT)) + 1;
    mtmp = mean(RR(indx,:));
    mE = mE + mtmp;
    sE = sE + mtmp.^2;
  end
  E = sqrt((sE/Nboot - mE.^2/Nboot^2));
end

% if adaptive warp sig so that on average the number of spikes
% under the kernel is the same but regions where there is 
% more data have a smaller kernel

if adapt 
  sigt = mean(R)*sig./R;
  RR = zeros(NT,N_pts);
  f = 1./(2*sigt.^2);
  for n=1:NT
    for m=1:ND(n)
      RR(n,:) = RR(n,:) + exp(-f.*(t-D(n,m)).^2);
    end
    RR(n,:) = RR(n,:).*(1./sqrt(2*pi*sigt.^2));
  end
  if NT > 1; R = mean(RR); else R = RR;end

  if err == 1
    E = sqrt(R./(2*NT*sigt*sqrt(pi)));
  elseif err == 2
    Nboot = 10;
    mE = 0;
    sE = 0;
    for b=1:Nboot
      indx = floor(NT*rand(1,NT)) + 1;
      mtmp = mean(RR(indx,:));
      mE = mE + mtmp;
      sE = sE + mtmp.^2;
    end
    E = sqrt((sE/Nboot - mE.^2/Nboot^2)); 
  end
end  

if plt == 'n';return;end
plot(t,R,plt)
hold on
if err > 0
  plot(t,R+2*E,'g')
  plot(t,R-2*E,'g')
end
%axis([T(1)+(4*sig) T(2)-(4*sig) 0 1.25*max(R)])
axis([T(1)+(4*sig) T(2)-(4*sig) 0 max(R+2*E)+10])
xlabel('time (s)')
ylabel('rate (Hz)')
title(['Trial averaged rate : Gaussian Kernel :'  ...
	    ' sigma = ' num2str(1000*sig) 'ms'])
hold off











