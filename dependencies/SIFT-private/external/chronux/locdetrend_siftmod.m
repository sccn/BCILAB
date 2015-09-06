function [data fitlines] =locdetrend(data,Fs,movingwin,method)
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
%  method       {'linear', 'constant'}. remove linear trend or mean
% Output:
% data:         (locally detrended data)
%
%
% NOTICE: This function is adapted -- for inclusion in SIFT -- from the open-source Chronux
% toolbox which can be downloaded from http://www.chronux.org/
%
% Modified by: Tim Mullen, 2011, SCCN/INC/UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

data = change_row_to_column(data);
[N,C]=size(data);
if nargin < 2 || isempty(Fs); Fs=1; end;
if nargin < 3 || isempty(movingwin); movingwin=[N/Fs N/Fs]; end;
if nargin < 4, method = 'linear'; end
Tw=movingwin(1); Ts=movingwin(2);
if Ts>Tw; error('Use step size shorter than window size'); end;
n=round(Fs*Tw);
dn=round(Fs*Ts);
if nargout>1
    fitlines = zeros(size(data));
end
if ~isreal(data) 
  yr=real(data); 
  yi=imag(data);
  if n==N;
     yr=detrend(yr,method);
     yi=detrend(yi,method);
     data=yr+1i*yi;
  else
     for ch=1:C
         
         tmp=runline_siftmod(yr(:,ch),n,dn,method); 
         yr=yr-tmp;
         tmp=runline_siftmod(yi(:,ch),n,dn,method); 
         yi=yi-tmp;
         data(:,ch)=yr+1i*yi;
     end;
  end;
else
  if n==N;
      if nargout>1
        datadet=detrend(data,method);
        fitlines = data-datadet;
        data = datadet;
      else
          data=detrend(data,method);
      end
  else
     for ch=1:C;
         tmp=runline_siftmod(data(:,ch),n,dn,method); 
         if nargout>1
             fitlines(:,ch) = tmp;
         end
         
         data(:,ch)=data(:,ch)-tmp;
     end;
  end
end
