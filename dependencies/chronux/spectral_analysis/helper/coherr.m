function [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave,numsp1,numsp2)
% Function to compute lower and upper confidence intervals on the coherency 
% given the tapered fourier transforms, errchk, trialave.
%
% Usage: [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave,numsp1,numsp2)
% Inputs:
% C     - coherence
% J1,J2 - tapered fourier transforms 
% err - [errtype p] (errtype=1 - asymptotic estimates; errchk=2 - Jackknife estimates; 
%                   p - p value for error estimates)
% trialave - 0: no averaging over trials/channels
%            1 : perform trial averaging
% numsp1    - number of spikes for data1. supply only if finite size corrections are required
% numsp2    - number of spikes for data2. supply only if finite size corrections are required
%
% Outputs: 
%          confC - confidence level for C - only for err(1)>=1
%          phistd - theoretical or jackknife standard deviation for phi for err(1)=1 and err(1)=2 
%                   respectively. returns zero if coherence is 1
%          Cerr  - Jacknife error bars for C  - only for err(1)=2

if nargin < 5; error('Need at least 5 input arguments'); end;
if err(1)==0; error('Need err=[1 p] or [2 p] for error bar calculation'); end;
if nargout==4  && err(1)==1; error('Cerr contains Jackknife errors: only computed when err(1) is 2'); end;
[nf,K,Ch]=size(J1);
errchk=err(1);
p=err(2);
pp=1-p/2;
%
% Find the number of degrees of freedom
%
if trialave;
   dim=K*Ch;
   dof=2*dim;
   dof1=dof;
   dof2=dof;
   Ch=1;
   if nargin>=6 && ~isempty(numsp1) 
      totspikes1=sum(numsp1);
      dof1=fix(2*totspikes1*dof/(2*totspikes1+dof));
   end
   if nargin==7 && ~isempty(numsp2); 
      totspikes2=sum(numsp2);
      dof2=fix(2*totspikes2*dof/(2*totspikes2+dof));
   end;
   dof=min(dof1,dof2);
   J1=reshape(J1,nf,dim);
   J2=reshape(J2,nf,dim);
else
   dim=K;
   dof=2*dim;
   dof1=dof;
   dof2=dof;
   for ch=1:Ch;
      if nargin>=6 && ~isempty(numsp1);
         totspikes1=numsp1(ch); 
        dof1=fix(2*totspikes1*dof/(2*totspikes1+dof));
      end;
      if nargin==7 && ~isempty(numsp2);
         totspikes2=numsp2(ch);
        dof2=fix(2*totspikes2*dof/(2*totspikes2+dof));
      end;
      dof(ch)=min(dof1,dof2);
   end;
end;
%
% variance of the phase
%
%
% Old code is the next few lines - new code is in the if statement below
% beginning line 87
%
% if isempty(find((C-1).^2 < 10^-16));
%    phierr = sqrt((2./dof(ones(nf,1),:)).*(1./(C.^2) - 1));  
% else
%    phierr = zeros(nf,Ch);
% end  

%
% theoretical, asymptotic confidence level
%
if dof <= 2
   confC = 1;
else     
   df = 1./((dof/2)-1);
   confC = sqrt(1 - p.^df);
end;
%
% Phase standard deviation (theoretical and jackknife) and jackknife
% confidence intervals for C
%
if errchk==1;
   totnum=nf*Ch;
   phistd=zeros(totnum,1); 
   CC=reshape(C,[totnum,1]); 
   indx=find(abs(CC-1)>=1.e-16);
   dof=repmat(dof,[nf,1]);
   dof=reshape(dof,[totnum 1]); 
   phistd(indx)= sqrt((2./dof(indx).*(1./(C(indx).^2) - 1))); 
   phistd=reshape(phistd,[nf Ch]);
elseif errchk==2;
    tcrit=tinv(pp,dof-1);
    for k=1:dim;
        indxk=setdiff(1:dim,k);
        J1k=J1(:,indxk,:);
        J2k=J2(:,indxk,:);
        eJ1k=squeeze(sum(J1k.*conj(J1k),2));
        eJ2k=squeeze(sum(J2k.*conj(J2k),2));
        eJ12k=squeeze(sum(conj(J1k).*J2k,2)); 
        Cxyk=eJ12k./sqrt(eJ1k.*eJ2k);
        absCxyk=abs(Cxyk);
        atanhCxyk(k,:,:)=sqrt(2*dim-2)*atanh(absCxyk);
        phasefactorxyk(k,:,:)=Cxyk./absCxyk;
%         indxk=setdiff(1:dim,k);
%         J1jk=J1(:,indxk,:);
%         J2jk=J2(:,indxk,:);
%         eJ1jk=squeeze(sum(J1jk.*conj(J1jk),2));
%         eJ2jk=squeeze(sum(J2jk.*conj(J2jk),2));
%         eJ12jk=squeeze(sum(conj(J1jk).*J2jk,2)); 
%         atanhCxyjk(k,:,:)=sqrt(2*dim-2)*atanh(abs(eJ12jk)./sqrt(eJ1jk.*eJ2jk));
    end; 
    atanhC=sqrt(2*dim-2)*atanh(C);
    sigma12=sqrt(dim-1)*squeeze(std(atanhCxyk,1,1));
%     sigma12=sqrt(dim-1)*squeeze(std(atanhCxyjk,1,1));
    if Ch==1; sigma12=sigma12'; end;
    Cu=atanhC+tcrit(ones(nf,1),:).*sigma12;
    Cl=atanhC-tcrit(ones(nf,1),:).*sigma12;
    Cerr(1,:,:) = tanh(Cl/sqrt(2*dim-2));
    Cerr(2,:,:) = tanh(Cu/sqrt(2*dim-2));
    phistd=(2*dim-2)*(1-abs(squeeze(mean(phasefactorxyk))));
    if trialave; phistd=phistd'; end;
end
% ncrit=norminv(pp);
% phierr=zeros([2 size(phistd)]); 
% phierr(1,:,:)=phi-ncrit*phistd; phierr(2,:,:)=phi+ncrit*phistd;
