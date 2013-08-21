function y_line=runline(y,n,dn)
% Running line fit (local linear regression)
%
% Usage: y_line=runline(y,n,dn);
%
% Inputs: 
% y: input 1-d time series (real)
% n: length of running window in samples
% dn: stepsize of window in samples
% 
% Outputs:
% y_line: local line fit to data
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
	y1=mean(tseg); 
	y2=mean((1:n)'.*tseg)*2/(n+1);
	a=(y2-y1)*6/(n-1); b=y1-a*(n+1)/2;
	yfit(j,:)=(1:n)*a+b;
	y_line((j-1)*dn+(1:n))=y_line((j-1)*dn+(1:n))+(yfit(j,:).*wt)';
	norm((j-1)*dn+(1:n))=norm((j-1)*dn+(1:n))+wt';
end
mask=find(norm>0); y_line(mask)=y_line(mask)./norm(mask);
indx=(nwin-1)*dn+n-1;
npts=length(y)-indx+1;
y_line(indx:end)=(n+1:n+npts)'*a+b;
