function [mask,Xdiff]=plotsigdiff(X1,X1err,X2,X2err,plt,t,f)
% Function to plot significant differences between two time-frequency arrays X1 and X2
% given errors X1err, X2err. 
% Usage: mask=plotsigdiff(X1,X1err,X2,X2err,plt,t,f)
%
% X1 err and X2err contain upper and lower confidence intervals for X1 and X2
% The plot generated is shows X1-X2 where the difference is significant
% either in dB or on a linear scale.
%
% Inputs:
% X1: input array t x f. Can also be a function of just the frequency. 
% X1err: lower and upper confidence intervals for X1: lower/upper x t x f
% X2: input array t x f. if vector then as row vector
% X2err: lower and upper condidence intervals for X2: lower/upper x t x f
% plt: 'l' for log, 'nl' for no log,'n' for no plot at all.
% t: t axis grid for plot. If X1,X2 are vectors, then specify t=1.
% f: f axis grid for plot.
%
% Outputs:
% mask: +1 for all t-f (or f) indices for which the X1 significantly greater than
% X2, -1 for all t-f (or f) indices for which X1 is significantly less than X2,
% and zero otherwise
%
% Xdiff: X1-X2
%
if nargin < 7; error('Need all arguments'); end;
% [T1,F1]=size(X1); [T2,F2]=size(X2); 
[T,F]=check_consistency(X1,X2);
if F==1;
    X1=X1'; X2=X2';F=length(X1); T=1;
end;
ystr='';
if T==1,
    mask=zeros(1,F);
    indxneg=find(X1<X2err(1,:) & X2>X1err(2,:));
    indxpos=find(X1>X2err(2,:) & X2<X1err(1,:));
    mask(indxneg)=-1;
    mask(indxpos)=+1;
    if strcmp(plt,'l'); 
        X1=10*log10(X1); X2=10*log10(X2); X1err=10*log10(X1err); X2err=10*log10(X2err);
        ystr= '  dB';
    end;
    subplot(311); plot(f,X1,f,X1err(1,:),f,X1err(2,:));
    title('Spectrum 1');
    xlabel('f')
    ylabel(['S1' ystr]);
    subplot(312); plot(f,X2,f,X2err(1,:),f,X2err(2,:));
    title('Spectrum 2');
    xlabel('f')
    ylabel(['S2' ystr]);
    subplot(313); plot(f,mask.*(X1-X2));
    title('Difference where significant');
    xlabel('f')
    ylabel(['S1-S2' ystr]);
else
    mask=zeros(T,F);
    for n=1:length(t);
        for m=1:length(f);
           if X1(n,m)<X2err(1,n,m) && X2(n,m)>X1err(2,n,m);
              mask(n,m)=-1;
           elseif X2(n,m)<X1err(1,n,m) && X1(n,m)>X2err(2,n,m);
              mask(n,m)=+1;
           end;
        end;
    end;
    if strcmp(plt,'l');
       X1=10*log10(X1);X2=10*log10(X2); %X1err=10*log10(X1err); X2err=10*log10(X2err);
        ystr='  dB';
    end;
    if ~strcmp(plt,'n');
        subplot(311); imagesc(t,f,X1'); axis xy; colorbar;
        xlabel('f')
        ylabel(['S1' ystr]);
        subplot(312); imagesc(t,f,X2'); axis xy; colorbar;
        xlabel('f')
        ylabel(['S2' ystr]);
	%     subplot(313); imagesc(t,f,(mask.*(X1-X2))'); axis xy; colorbar
        subplot(313); imagesc(t,f,mask'); axis xy; colorbar
        xlabel('f')
        ylabel('Significance');
    end;
end
Xdiff=X1-X2;
