function plot_matrix(X,t,f,plt,Xerr)
% Function to plot a time-frequency matrix X. Time and frequency axes are in t and f.
% If error bars are specified in Xerr,
% it also plots them. Xerr contains upper and lower confidence intervals 
% on X. 
% Usage: plot_matrix(X,t,f,plt,Xerr)
% Inputs:
% X: input vector as a function of time and frequency (t x f)
% t: t axis grid for plot. Default [1:size(X,1)]
% f: f axis grid for plot. Default. [1:size(X,2)]
% plt: 'l' for log, 'n' for no log.
% Xerr: lower and upper confidence intervals for X1: lower/upper x t x f.
if nargin < 1; error('Need data'); end;
[NT,NF]=size(X);
if nargin < 2;
    t=1:NT;
end;
if nargin < 3;
    f=1:NF;
end;
if length(f)~=NF || length(t)~=NT; error('axes grid and data have incompatible lengths'); end;
if nargin < 4 || isempty(plt);
    plt='l';
end;
if strcmp(plt,'l');
    X=10*log10(X);
    if nargin ==5; Xerr=10*log10(Xerr); end;
end;

if nargin < 5;
   imagesc(t,f,X');axis xy; colorbar; title('Spectrogram');
else
   subplot(311); imagesc(t,f,squeeze(Xerr(1,:,:))'); axis xy; colorbar; title('Lower confidence');
   subplot(312); imagesc(t,f,X'); title('X');axis xy; colorbar;
   subplot(313); imagesc(t,f,squeeze(Xerr(2,:,:))'); axis xy; colorbar; title('Upper confidence');
end;
xlabel('t');ylabel('f');

