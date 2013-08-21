function  [C,phi,S12,confC,phierr,Cerr]=cohmathelper(J,err,Nsp)
% Helper function called by coherency matrix computations.
%
% Usage: [C,phi,S12,confC,phierr,Cerr]=cohmathelper(J,err,Nsp)
% Inputs:
% J   : Fourier transforms of data
% err : [0 p] or 0 for no errors; [1 p] for theoretical confidence level, 
%       [2 p] for Jackknife (p - p value)
% Nsp : pass the number of spikes in each channel if finite size corrections are desired
%
% Outputs:
%
% C   : coherence
% phi : phase of coherency
% S12 : cross spectral matrix
% confC : confidence level for coherency - only for err(1)>=1
%       phierr - standard deviation for phi (note that the routine gives phierr as phierr(1,...) 
%                and phierr(2,...) in order to incorporate Jackknife (eventually). 
%                Currently phierr(1,...)=phierr(2,...). Note that phi + 2 phierr(1,...) and phi -2 
%                phierr(2,...) will give 95% confidence bands for phi - only for err(1)>=1
% Cerr  : error bars for coherency (only for Jackknife estimates)-only for err(1)=2
%

errtype=err(1);
trialave=0;
[nf,K,Ch]=size(J);
clear K
confC=zeros(Ch,Ch);
C=zeros(nf,Ch,Ch);
S12=zeros(nf,Ch,Ch);
phi=zeros(nf,Ch,Ch);
phierr=zeros(2,nf,Ch,Ch);
if errtype==2; Cerr=zeros(2,nf,Ch,Ch);end;

for ch1=1:Ch;
     J1=squeeze(J(:,:,ch1));
     C(1:nf,ch1,ch1)=1;
     phi(1:nf,ch1,ch1)=0;
%      if errtype==2; 
%           phierr(1:nf,ch1,ch1)=0;
%           Cerr(1:2,1:nf,ch1,ch1)=0;
%      elseif errtype==1
%            phierr(1:2,1:nf,ch1,ch1)=0;
%      end;
     s1=squeeze(mean(conj(J1).*J1,2));
     for ch2=1:ch1-1;
          J2=squeeze(J(:,:,ch2));
          s12=squeeze(mean(conj(J1).*J2,2));
          s2=squeeze(mean(conj(J2).*J2,2));
          C12=s12./sqrt(s1.*s2);
          C(:,ch1,ch2)=abs(C12);
          C(:,ch2,ch1)=C(:,ch1,ch2);
          phi(:,ch1,ch2)=angle(C12);
          phi(:,ch2,ch1)=phi(:,ch1,ch2);
          S12(:,ch1,ch2)=s12;
          S12(:,ch2,ch1)=S12(:,ch1,ch2);
          if errtype==2 
             if nargin<3;
                 [conf,phie,Ce]=coherr(abs(C12),J1,J2,err,trialave);
             else
                 [conf,phie,Ce]=coherr(abs(C12),J1,J2,err,trialave,Nsp(ch1),Nsp(ch2));
             end
             confC(ch1,ch2)=conf; 
             phierr(1,:,ch1,ch2)=phie;phierr(2,:,ch1,ch2)=phie;
             Cerr(1,:,ch1,ch2)=Ce(1,:);
             Cerr(2,:,ch1,ch2)=Ce(2,:);
             confC(ch2,ch1)=conf; 
             phierr(1,:,ch2,ch1)=phie;phierr(2,:,ch2,ch1)=phie;
             Cerr(:,:,ch2,ch1)=Ce;
          elseif errtype==1
             if nargin<3;
                 [conf,phie]=coherr(abs(C12),J1,J2,err,trialave);
             else
                 [conf,phie]=coherr(abs(C12),J1,J2,err,trialave,Nsp(ch1),Nsp(ch2));
             end
             confC(ch1,ch2)=conf; 
             phierr(1,:,ch1,ch2)=phie;phierr(2,:,ch1,ch2)=phie;
             confC(ch2,ch1)=conf; 
             phierr(1,:,ch2,ch1)=phie;phierr(2,:,ch2,ch1)=phie;
          end;
     end;
end;
