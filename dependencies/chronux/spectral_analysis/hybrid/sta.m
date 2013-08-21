function[s,t,E] = sta(data_spk,data_lfp,smp,plt,w,T,D,err)
% Spike Triggered Average                            
%     Usage: [s,t,E] = sta(data_spk,data_lfp,smp,plt,w,T,D,err)
%     
% Inputs                                              
%                                                    
% Note that all times have to be consistent. If data_spk
% is in seconds, so must be sig and t. If data_spk is in 
% samples, so must sig and t. The default is seconds.
%
% data_spk    - strucuture array of spike times data 
%               or NaN padded matrix
% data_lfp    - array of lfp data(samples x trials)         
%                                                    
% Optional...                                        
% plt 'n'|'r' etc                                    
% width kernel smoothing in s                        
% T = [-0.1 0.1] - extract this range about each spk 
% D = plot spike triggered average out to [D1 D2]    
% err = calcluate error bars (bootstrap)             
%                                                    
% Outputs:                                             
%                                                    
% s  spike triggered average                         
% t  times                                           
% E  bootstrap standard err                          

if nargin < 3;error('Require spike, lfp and lfp times');end
if isstruct(data_spk)
   [data_spk]=padNaN(data_spk); % create a zero padded data matrix from input structural array
   sz=size(data_spk); 
   if sz(1)>sz(2); data_spk=data_spk'; end;% transpose data to get in form compatible with Murray's routine
else
   sz=size(data_spk);
   if sz(1)>sz(2); data_spk=data_spk'; end;% transpose data to get in form compatible with Murray's routine
end
sz=size(data_lfp);
if sz(1)>sz(2); data_lfp=data_lfp'; end;% transpose data to get in form compatible with Murray's routine
verbose = 1;
t = smp;
if nargin < 4; plt = 'r'; end
if nargin < 5; w = 0; end
if nargin < 6; T = [min(t) max(t)]; end
if nargin < 7; D = 0.25*[-1 1]; end
if nargin < 8; err = 1;end

if isempty(plt); plt = 'r'; end
if isempty(w); w = 0; end
if isempty(T); T = [min(t) max(t)]; end
if isempty(D); D = 0.25*[-1 1]; end
if isempty(err); err = 1;end

if w > (T(2)-T(1))/2
  disp('Smoothing > data segment : should be in seconds : turn off smoothing')
  w = 0;
end

sz = size(data_spk);
NT = sz(1);
mlfp = 0;
slfp = 0;
Nspk = 0;
smp = t(2)-t(1);
if D(1) <= 0 && D(2) >= 0
  t1 = [D(1):smp:(-smp+eps) 0:smp:D(2)+eps];
else
  t1 = (round(D(1)/smp)*smp):smp:D(2);
end

% count up the spikes...

if err
  for n=1:NT
    indx = find(t>T(1)&t<T(2));
%     lfp = data_lfp(n,indx);
    spk = data_spk(n,data_spk(n,:)>T(1) && data_spk(n,:)<T(2) && data_spk(n,:)~=0);
    tt = t(indx);
    if ~isempty(spk) > 0
      ND = length(spk);
      for s=1:ND
        spktime = spk(s);      
        t0 = tt-spktime + eps;
        if min(t0)>(D(1)-smp) || max(t0)<(D(2)+smp); break;end
        Nspk = Nspk + 1;
      end
    end
  end
  Err = zeros(Nspk,length(t1));
  Nspk = 0;
end

for n=1:NT
  indx = find(t>T(1)&t<T(2));
  lfp = data_lfp(n,indx);
  spk = data_spk(n,data_spk(n,:)>T(1) && data_spk(n,:)<T(2) && data_spk(n,:)~=0);
  tt = t(indx);
  if ~isempty(spk) > 0
    ND = length(spk);
    for s=1:ND
      spktime = spk(s);      
      t0 = tt-spktime + eps;
      if min(t0) < (D(1)-smp) && max(t0) > (D(2)+smp);
        indx = find(t0<D(1));
        indx = indx(length(indx));   
        offset = (t0(indx)-D(1))/smp; 
        indx = indx:(indx+length(t1)-1);
        lfp_t1 = lfp(indx) + (lfp(indx+1)-lfp(indx))*offset;
        Nspk = Nspk + 1;
        mlfp = mlfp + lfp_t1;
        slfp = slfp + lfp_t1.^2;      
        Err(Nspk,:) = lfp_t1;
      end
    end    
  end  
end
if Nspk == 0
  if verbose;disp('No spikes in interval');end
  t = t1;
  s = zeros(length(t),1);
  E = zeros(length(t),1);
  return
end
mlfp = mlfp/Nspk;
slfp = slfp/Nspk;
stdlfp = sqrt((slfp - mlfp.^2)/Nspk);

% local smoother...

N = fix(w/smp);
if N > 5
  mlfp = locsmooth(mlfp,N,fix(N/2)); 
end

% bootstrap errorbars...

if err == 1;
  Nboot = 20;  
  blfp = 0;
  slfp = 0;
  for n = 1:Nboot
    indx = floor(Nspk*rand(1,Nspk)) + 1;
    lfptmp = mean(Err(indx,:));
    if N > 5
      lfptmp = locsmooth(lfptmp,N,fix(N/2));
    end
    blfp = blfp + lfptmp;
    slfp = slfp + lfptmp.^2;
  end
  stdlfp = sqrt((slfp/Nboot - blfp.^2/Nboot^2));
end  

s = mlfp-mean(mlfp);
E = stdlfp;
t = t1;

%cols = 'krbgycm';
if plt == 'n';return;end
plot(t1,s,plt)
xax = get(gca,'xlim');
%yax = get(gca,'ylim');
if err == 1
  me = real(2*mean(stdlfp));
  line(xax,me*[1 1],'color','b')
  line(xax,-me*[1 1],'color','b')
  line(xax,0*[1 1],'color','k')
%  errorbar(0.1*xax(2)+0.9*xax(1),0.1*yax(2)+0.9*yax(1), ...
%         mean(stdlfp),'k')
%plot(0.1*xax(2)+0.9*xax(1),0.1*yax(2)+0.9*yax(1),'k.')
end
     
title(['spike triggered average : ' ...
      num2str(Nspk) ' used : '     ...
      ' Errorbars are two std err']);
%line(get(gca,'xlim'),mean(mlfp)*[1 1],'color','k')  
line([0 0],get(gca,'ylim'),'color','k')
hold off








