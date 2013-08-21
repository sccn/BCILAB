function [mintime, maxtime]=minmaxsptimes(data)
% Find the minimum and maximum of the spike times in each channel
% Usage: [mintime, maxtime]=minmaxsptimes(data)
% Input:
% data  (spike times as a structural array of multiple dimensions e.g. channels; channels x trials; 
%        can also accept a 1d matrix of spike times)
% Output:
% mintime       (minimum of the spike time across channels)
% maxtime       (maximum of the spike time across channels)
%
dtmp='';
if isstruct(data)
   data=reshape(data,numel(data),1);
   C=size(data,1);
   fnames=fieldnames(data);
   mintime=zeros(1,C); maxtime=zeros(1,C);
   for ch=1:C
     eval(['dtmp=data(ch).' fnames{1} ';'])
     if ~isempty(dtmp)
        maxtime(ch)=max(dtmp);
        mintime(ch)=min(dtmp);
     else
        mintime(ch)=NaN;
        maxtime(ch)=NaN;
     end
   end;
   maxtime=max(maxtime); % maximum time
   mintime=min(mintime); % minimum time
else
     dtmp=data;
     if ~isempty(dtmp)
        maxtime=max(dtmp);
        mintime=min(dtmp);
     else
        mintime=NaN;
        maxtime=NaN;
     end
end
if mintime < 0 
   error('Minimum spike time is negative'); 
end
