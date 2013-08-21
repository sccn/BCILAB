function plotsig(C,sig,t,f,c)
% Function to plot C where it is higher than a threshold sig
% useful for plotting coherence
% Usage: plotsig(C,sig,t,f)
% Inputs:
% C: input array t x f - also works for a single vector
% sig: significance level
% t: t axis grid for plot
% f: f axis grid for plot.
% c: color to use (default blue)-only meaningful for a line plot
if nargin < 4; error('Need at least 4 arguments'); end;
if nargin < 5 | isempty(c); c='b'; end;
[T,F]=size(C);
if F==1; C=C'; [T,F]=size(C);end;
if T~=length(t) | F~=length(f);
    error('frequency and/or time axes are incompatible with data'); 
end;
if T==1;
    dim=max(T,F);
    C=C(:);
    indx=find(C>sig);
    plot(f,C,c); 
%     plot(f,C,f,mask.*C)
    line(get(gca,'xlim'),[sig sig]);
    xlabel('f'); ylabel('|C|');
else
    mask=zeros(T,F);
    for n=1:length(t);
        for m=1:length(f);
           if C(n,m)>sig
              mask(n,m)=1;
           end;
        end;
    end;
    imagesc(t,f,(mask.*C)'); axis xy; colorbar
    xlabel('t'); ylabel('f');
end;