function data=createdatamatpb(data,E,Fs,win)
%
% Helper function to create an event triggered matrix from a single
% channel of data. 
% Usage: data=createdatamatpb(data,E,Fs,win)
% Inputs:
% data   (input time series as a single vector) - required
% E      (events to use as triggers) - required 
% Fs     (sampling frequency of data) - required
% win    (window around triggers to use data matrix -[winl winr]) - required 
%          e.g [1 1] uses a window starting 1 sec before E and
%              ending 1 sec after E if E is in secs
% Note that E, Fs, and win must have consistent units 
%
% Outputs:
% data      (event triggered data)
%
if nargin < 4
   error('Need all input arguments'); 
end
NE=length(E);
nwinl=round(win(1)*Fs);
nwinr=round(win(2)*Fs);
nE=floor(E*Fs)+1;
datatmp=[];
for n=1:NE;
    indx=nE(n)-nwinl:nE(n)+nwinr-1; 
    datatmp=[datatmp data(indx)];
end
data=datatmp;