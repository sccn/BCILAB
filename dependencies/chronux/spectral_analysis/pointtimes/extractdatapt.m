function data=extractdatapt(data,t,offset)
% Extract segements of spike times between t(1) and t(2)
% Usage: data=extractdatapt(data,t,offset)
%
% Input:
% data: structural array of spike times for each channel/trial or a single
%       array of spike times
% t   : time as a 2d vector [t(1) t(2)]
% offset: 0/1 - if 1, store the spike times relative to start of window i.e. t(1)
%         if 0, don't reset the times. Default 0. 
% Note that all times can be in arbitrary units. But the units have to be
% consistent. So, if E is in secs, win, t have to be in secs, and Fs has to
% be Hz. If E is in samples, so are win and t, and Fs=1. In case of spike
% times, the units have to be consistent with the units of data as well.
%
% Output:
% data: spike times between t(1) and t(2)
if nargin < 2; error('Need data and times'); end;
if t(1) < 0 || t(2)<=t(1);
    error('times cannot be negative and t(2) has to greater than t(1)');
end;
if nargin < 3 || isempty(offset); offset=0; end;
if isstruct(data); 
    C=length(data);
elseif min(size(data))~=1; 
    error('Can only accept single vector data unless it is a struct array'); 
else
    C=1;
    data=change_row_to_column(data);
end;
%fnames=fieldnames(data);
d2(1:C)=struct('times',[]);
for c=1:C,
    if isstruct(data)
       fnames=fieldnames(data);
       eval(['dtmp=data(c).' fnames{1} ';'])
    else
       dtmp=data(:);
    end
%     eval(['dtmp=data(c).' fnames{1} ';' ])
    sp=dtmp(dtmp>=t(1) & dtmp<t(2));
    if offset==1; d2(c).times=sp-t(1); 
    else d2(c).times=sp;end
end;
data=d2;
