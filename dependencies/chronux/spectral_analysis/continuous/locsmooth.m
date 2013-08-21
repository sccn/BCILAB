function data=locsmooth(data,Fs,Tw,Ts)
%  Running line fit (using local linear regression) - 1d only, continuous
%  processes
%  Usage: data=locsmooth(data,Fs,Tw,Ts)
%  Inputs:
% Note that units of Fs, movinwin have to be consistent.
%  data  (single vector) 
%  Fs    (sampling frequency) - optional. Default 1
%  Tw    (length of moving window) - optional.  Default. full length of data (global detrend)
%  Ts    (step size) - optional. Default Tw/2.
% 
% Output:
% data   (locally smoothed data).
data=change_row_to_column(data);
N=size(data,1);
if nargin < 2; Fs=1; end;
if nargin < 3; Tw=N/Fs; end;
if nargin < 4; Ts=Tw/2; end;

n=round(Fs*Tw);
dn=round(Fs*Ts);
if ~isreal(data) 
  yr=real(data); 
  yi=imag(data); 
  tmp=runline(yr,n,dn); 
  yr=tmp';
  tmp=runline(yi,n,dn); 
  yi=tmp';
  data=yr+i*yi;
else
  tmp=runline(data,n,dn); 
  data=tmp';
end
