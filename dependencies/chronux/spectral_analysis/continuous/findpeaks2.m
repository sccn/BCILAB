function xmax=findpeaks(data,threshold)
% Helper function to find peaks in a given continuous valued time series x
% Usage: xmax=findpeaks(data,threshold)
% Input:
%      data     (data in time x channels/trials form or a single vector)
%      threshold (if specified returns locations of peaks at which data exceeds threshold) - optional
% Output:
%      xmax     (locations of local maxima of data in a structure array of dimensions channels/trials)
if nargin < 1; error('Need data'); end;
data=change_row_to_column(data);
C=size(data,2);
pp1=[data(1,:);data(1:end-1,:)];
pp2=[data(2:end,:);data(end,:)];
xmax(1:C)=struct('loc',[]);
for ch=1:C,
   if nargin ==1
     xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>0 & data(:,ch)-pp2(:,ch)>0)];
   else
     xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>0 & data(:,ch)-pp2(:,ch)>0 & data(:,ch)>threshold)];
   end
end