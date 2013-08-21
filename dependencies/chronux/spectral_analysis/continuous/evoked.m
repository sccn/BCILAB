function [V,t,Err] = evoked(data,Fs,win,width,plt,err)
% Function to calculate the evoked response given continuous data in the
% form time x channels
% Usage [V,t,Err] = evoked(data,Fs,win,width,plt,err)
% 
% Inputs  
%   Note that all times can be in arbitrary units. But the units have to be
%   consistent. So, if win is in secs, width is in secs and Fs has to be Hz. 
%   If win is in samples, so is width and Fs=1.
%
%    data(times, channels/trials or a single vector)      (required)    
%    Fs  sampling frequency            (required)
%    win   subsection of data to be used. Default all available data
%    width (s) of smoothing kernel. Default 50 samples                  
%    plt plot 'n' for no plot, otherwise plot color. Default blue colored lines.                                  
%    err = 0/1. Default 1=calculate bootstrap errorbars.                     
%                                                         
% Outputs                                             
%    V = evoked potential                                 
%    t = times of evaluation                              
%    Err = bootstrap statdard deviation                   

if nargin < 2;error('Data, sampling frequency required');end
data=change_row_to_column(data);
N=size(data,1);
data=data';
if nargin <3; win = [0 (N-1)/Fs];end
if nargin <4; width = 50/Fs;end
if nargin <5; plt = 'b';end
if nargin <6;err = 1;end
T=win;
if isempty(T); T = [0 (N-1)/Fs];end
if isempty(width); width = 50/Fs;end
if isempty(plt); plt = 'b';end
if isempty(err);err = 1;end

t = min(T):1/Fs:max(T);
if nargin >= 5
  indx = find(t>T(1) & t<T(2));
  t = t(indx);
  data = data(:,indx);
end

if width > (t(length(t))-t(1))/2
  disp('Width is too large for data segment: should be in seconds')
  disp('Turn off smoothing')
  width = 0;
end

s = t(2)-t(1);
N = fix(width/s);
NT = length(data(:,1));

if NT > 1
    mdata = mean(data);
else
    mdata = data;
end
if N > 4
  smdata = locsmooth(mdata,N,fix(N/2)); 
else
  smdata = mdata;  
end
  
% if errorbars requested then do a bootstrap over trials...

Err = 0;
if NT < 4; 
  disp('Too few trials: no errorbars calculated')
  err = 0;    
end

if err ~= 0 && NT > 1
  Nboot = 10;
  bevk = 0;
  sevk = 0;
  for b=1:Nboot
    indx = floor(NT*rand(1,NT)) + 1;
    evktmp = mean(data(indx,:));
    if N > 4
      evktmp = locsmooth(evktmp,N,fix(N/2));
    end
    bevk = bevk + evktmp;
    sevk = sevk + evktmp.^2;
  end
  stdevk = sqrt((sevk/Nboot - bevk.^2/Nboot^2));
  Err = stdevk;
end

V = smdata;
if plt ~= 'n'
  plot(t,smdata,plt)
  hold on
  mn = mean(smdata);
  ax = get(gca,'xlim');
  line(ax,mn*[1 1],'color','k')
  if err
    line(ax,(mn+2*mean(stdevk))*[1 1],'color','r')
    line(ax,(mn-2*mean(stdevk))*[1 1],'color','r')
    hold off
  end
end









