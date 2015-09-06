function y_line=runline_siftmod(y,n,dn,method)
% Running line fit (local linear regression)
%
% Usage: y_line=runline(y,n,dn);
%
% Inputs: 
% y: input 1-d time series (real)
% n: length of running window in samples
% dn: stepsize of window in samples
% method: 'linear' or 'constant' (remove linear fit or mean)
%
% Outputs:
% y_line: local line fit to data
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

y=y(:);
nt=length(y);
y_line=zeros(nt,1);
norm=y_line;
nwin=ceil((nt-n)/dn);
yfit=zeros(nwin,n);
xwt=((1:n)-n/2)/(n/2);
wt=(1-abs(xwt).^3).^3;
for j=1:nwin, 
	tseg=y(dn*(j-1)+1:dn*(j-1)+n);
	y1=sum(tseg)/n;
    
    if strcmpi(method,'linear')
        y2=sum((1:n)'.*tseg)*2/(n*(n+1));
        a=(y2-y1)*6/(n-1); b=y1-a*(n+1)/2;
        yfit(j,:)=(1:n)*a+b;
    else
        yfit(j,:)=y1(ones(size(tseg)));
    end
    y_line((j-1)*dn+(1:n))=y_line((j-1)*dn+(1:n))+(yfit(j,:).*wt)';
	norm((j-1)*dn+(1:n))=norm((j-1)*dn+(1:n))+wt';
end
mask=find(norm>0); y_line(mask)=y_line(mask)./norm(mask);
indx=(nwin-1)*dn+n-1;
npts=length(y)-indx+1;
if strcmpi(method,'linear')
    y_line(indx:end)=(n+1:n+npts)'*a+b;
end
