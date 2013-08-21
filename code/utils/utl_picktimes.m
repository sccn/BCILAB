function Y = utl_picktimes(X,wnd)
% Average the input array within the given index ranges, for each row and plane.
% Result = utl_picktimes(Data, Windows)
%
% In:
%   Data    : [Channels x Samples x Trials] data array
%   Windows : [Ranges x 2] array specifying the beginning and end of each range, in samples
%
% Out:
%   Result  : averaged sub-ranges of the original data array; sized [Channels x Ranges x Trials]
%
% Examples:
%   % for a given 3d data array, average intervals of 20:30 samples, 50:100 samples, and 100:200 samples
%   % for each epoch and each channel
%   averages = utl_picktimes(EEG.data,[20 30;50 100;100 200])
%
% See also:
%   set_picktimes
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-20

[C,S,T] = size(X);

% fix the windows, if necessary
wnd = min(max(wnd,1),S);
wnd(:,2) = wnd(:,2)+1;
wnd = wnd(wnd(:,2)>wnd(:,1),:);
W = size(wnd,1);

% calc offset & coverage values for first & last sample
fo = floor(wnd(:,1));
fc = min(fo+1,wnd(:,2)) - wnd(:,1);
lo = ceil(wnd(:,2)-1);
lc = wnd(:,2) - max(lo,wnd(:,1));
% calc fully overlapped sample range, inverse length
for r=1:W
    full{r} = fo(r)+1:lo(r)-1; end
ilen = 1./(fc+cellfun(@length,full)'+lc); 

Y = zeros(C,W,T);
% accumulate for every range...
for r=1:W    
    Y(:,r,:) = (X(:,fo(r),:)*fc(r)+sum(X(:,full{r},:),2)+X(:,lo(r),:)*lc(r))*ilen(r); end
