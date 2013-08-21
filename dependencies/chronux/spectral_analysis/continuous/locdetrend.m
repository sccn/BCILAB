function data=locdetrend(data,Fs,movingwin)
%  Remove running line fit (using local linear regression)-continuous
%  processes
%  Usage: data=locdetrend(data,Fs,movingwin)
%  Inputs:
%  Note that units of Fs, movinwin have to be consistent.
%  data         (data as a matrix times x channels or a single vector) 
%  Fs           (sampling frequency) - optional. Default 1
%  movingwin    (length of moving window, and stepsize) [window winstep] - optional.
%                   Default. window=full length of data (global detrend).
%                   winstep=window -- global detrend
% 
% Output:
% data:         (locally detrended data)
data=change_row_to_column(data);
[N,C]=size(data);
if nargin < 2 || isempty(Fs); Fs=1; end;
if nargin < 3 || isempty(movingwin); movingwin=[N/Fs N/Fs]; end;
Tw=movingwin(1); Ts=movingwin(2);
if Ts>Tw; error('Use step size shorter than window size'); end;
n=round(Fs*Tw);
dn=round(Fs*Ts);
if ~isreal(data) 
  yr=real(data); 
  yi=imag(data);
  if n==N;
     yr=detrend(yr);
     yi=detrend(yi);
     data=yr+i*yi;
  else;
     for ch=1:C
         tmp=runline(yr(:,ch),n,dn); 
         yr=yr-tmp;
         tmp=runline(yi(:,ch),n,dn); 
         yi=yi-tmp;
         data(:,ch)=yr+i*yi;
     end;
  end;
else
  if n==N;
     data=detrend(data);
  else;
     for ch=1:C;
         tmp=runline(data(:,ch),n,dn); 
         data(:,ch)=data(:,ch)-tmp;
     end;
  end
end
