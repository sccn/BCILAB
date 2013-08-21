function[S,tau,tc] = staogram(data_spk,data_lfp,smp,plt,Tc,Tinc,Tw,w,D)
%
% staogram : calculates a moving window spike triggered ave %
%   Usage:[S,tau,tc] = staogram(data_spk,data_lfp,smp,plt,Tc,Tinc,Tw,w,D)
%   
%                 ******** INPUT *********                  
% Note that all times have to be consistent. If data_spk
% is in seconds, so must be sig and t. If data_spk is in 
% samples, so must sig and t. The default is seconds.
%
% data_spk    - strucuture array of spike times data        
% data_lfp    - array of lfp data(samples x trials)         
% smp         - lfp times of samples                        
%                                                           
% Optional...                                               
%                                                           
% Parameter                                                 
%                                                           
% plt     'y'|'n'                                           
%                                                           
% 'y' standard staogram                                     
% 'n' no plot                                               
%                                                           
% Tc = start and end times (centres)           whole trial       
% Tinc = time increment between windows             0.1            
% Tw = time window width                          0.3            
% w = smoothing width in seconds                            
% D = plot sta out to on axis [D(1) D(2)] s                 
%                                                           
%                ******** OUTPUT ********                   
%  S spike triggered average                                
%  tau - lag                                                
%  tc  - bin centers                                        


% setup defaults...
if nargin < 3;error('Require spike, lfp and lfptimes ');end
[data_spk]=padNaN(data_spk); % create a zero padded data matrix from input structural array
data_spk=data_spk'; % transposes data to get it in a form compatible with Murray's routine
if nargin < 4; plt = 'y';end
if nargin < 6; Tinc = 0.1; end
if nargin < 7; Tw = 0.5;end
if nargin < 8; w = 0.01;end
if nargin < 9; D = 0.15*[-1 1]; end
if nargin < 5; 
    Tc(1) = min(data_spk(:,1)) + Tw/2;
    Tc(2) = max(max(data_spk)) - Tw/2;
end

if isempty(plt); plt = 'y';end
if isempty(Tinc); Tinc = 0.1; end
if isempty(Tw); Tw = 0.5;end
if isempty(w); w = 0.01;end
if isempty(D); D = 0.15*[-1 1]; end
if isempty(Tc); 
    Tc(1) = min(data_spk(:,1)) + Tw/2;
    Tc(2) = max(max(data_spk)) - Tw/2;
end


%  round to nearest tinc...

t = smp;
Tc(1) = ceil(Tc(1)/Tinc)*Tinc;
Tc(2) = floor(Tc(2)/Tinc)*Tinc;
tc = Tc(1):Tinc:Tc(2);
for tt=1:length(tc)
  T = [tc(tt)-Tw/2 tc(tt)+Tw/2];
  if tt == 1
    [SS,tau] = sta(data_spk,data_lfp,t,'y',w,T,D,0);
    S = zeros(length(tc),length(SS));
  else
    [SS,tau] = sta(data_spk,data_lfp,t,'y',w,T,D,0);
  end
  S(tt,:) = SS';
  S(tt,:) = SS';
end

if ~strcmp(plt,'n')
  imagesc(tc,tau,squeeze(S)')
  set(gca,'ydir','normal')
  xlabel('time (s)')
  ylabel('frequency (Hz)')
  colorbar;
%  axes(h)
%  line(get(h,'xlim'),conf_C*[1 1],'color','k','linewidth',5)
end










