function [dz,vdz,Adz]=two_group_test_coherence(J1c1,J2c1,J1c2,J2c2,p,plt,f)
% function [dz,vdz,Adz]=two_group_test_coherence(J1c1,J2c1,J1c2,J2c2,p)
% Test the null hypothesis (H0) that data sets J1c1,J2c1,J1c2,J2c2 in 
% two conditions c1,c2 have equal population coherence
%
% Usage:
% [dz,vdz,Adz]=two_sample_test_coherence(J1c1,J2c1,J1c2,J2c2,p)
%
% Inputs:
% J1c1   tapered fourier transform of dataset 1 in condition 1
% J2c1   tapered fourier transform of dataset 1 in condition 1
% J1c2   tapered fourier transform of dataset 1 in condition 2
% J2c2   tapered fourier transform of dataset 1 in condition 2
% p      p value for test (default: 0.05)
% plt    'y' for plot and 'n' for no plot
% f      frequencies (useful for plotting)
%
%
% Dimensions: J1c1,J2c2: frequencies x number of samples in condition 1
%              J1c2,J2c2: frequencies x number of samples in condition 2
%              number of samples = number of trials x number of tapers
% Outputs:
% dz    test statistic (will be distributed as N(0,1) under H0
% vdz   Arvesen estimate of the variance of dz
% Adz   1/0 for accept/reject null hypothesis of equal population
%       coherences based dz ~ N(0,1)
% 
% Note: all outputs are functions of frequency
%
% References: Arvesen, Jackkknifing U-statistics, Annals of Mathematical
% Statisitics, vol 40, no. 6, pg 2076-2100 (1969)


if nargin < 4; error('Need four sets of Fourier transforms'); end;
if nargin < 6 || isempty(plt); plt='n'; end;

% 
% Test for matching dimensionalities
%
if size(J1c1)~=size(J2c1) | size(J1c2)~=size(J2c2) | size(J1c1,1)~=size(J1c2,1);
    error('Need matching dimensionalities for the Fourier transforms: Check the help file for correct dimensionalities');
else;
    m1=size(J1c1,2); % number of samples, condition 1
    m2=size(J1c2,2); % number of samples, condition 2
    dof1=2*m1; % number of degrees of freedom in the first condition estimates
    dof2=2*m2; % number of degrees of freedom in the second condition estimates
end;
if nargin < 7 || isempty(f); f=size(J1c1,1); end;
if nargin < 5; p=0.05; end; % set the default p value

%
% Compute the individual condition spectra, coherences
%
S12c1=conj(J1c1).*J2c1; % individual sample cross-spectrum, condition 1
S12c2=conj(J1c2).*J2c2; % individual sample cross-spectrum, condition 2
S1c1=conj(J1c1).*J1c1; % individual sample spectrum, data 1, condition 1
S2c1=conj(J2c1).*J2c1; % individual sample spectrum, data 2, condition 1
S1c2=conj(J1c2).*J1c2; % individual sample spectrum, data 1, condition 2
S2c2=conj(J2c2).*J2c2; % individual sample spectrum, data 2, condition 2

Sm12c1=squeeze(mean(S12c1,2)); % mean cross spectrum, condition 1
Sm12c2=squeeze(mean(S12c2,2)); % mean cross spectrum, condition 2
Sm1c1=squeeze(mean(S1c1,2)); % mean spectrum, data 1, condition 1
Sm2c1=squeeze(mean(S2c1,2)); % mean spectrum, data 2, condition 1
Sm1c2=squeeze(mean(S1c2,2)); % mean spectrum, data 1, condition 1
Sm2c2=squeeze(mean(S2c2,2)); % mean spectrum, data 2, condition 1

Cm12c1=abs(Sm12c1./sqrt(Sm1c1.*Sm2c1)); % mean coherence, condition 1
Cm12c2=abs(Sm12c2./sqrt(Sm1c2.*Sm2c2)); % mean coherence, condition 2

Ccm12c1=Cm12c1; % mean coherence saved for output
Ccm12c2=Cm12c2; % mean coherence saved for output
%
% Compute the statistic dz, and the probability of observing the value dz
% given an N(0,1) distribution i.e. under the null hypothesis
%
z1=atanh(Cm12c1)-1/(dof1-2); % Bias-corrected Fisher z, condition 1
z2=atanh(Cm12c2)-1/(dof2-2); % Bias-corrected Fisher z, condition 2
dz=(z1-z2)/sqrt(1/(dof1-2)+1/(dof2-2)); % z statistic
%
% The remaining portion of the program computes Jackknife estimates of the mean (mdz) and variance (vdz) of dz
% 
samples1=[1:m1];
samples2=[1:m2];
%
% Leave one out of one sample
%
for i=1:m1;
    ikeep=setdiff(samples1,i); % all samples except i
    Sm12c1=squeeze(mean(S12c1(:,ikeep),2)); % 1 drop mean cross-spectrum, condition 1
    Sm1c1=squeeze(mean(S1c1(:,ikeep),2)); % 1 drop mean spectrum, data 1, condition 1
    Sm2c1=squeeze(mean(S2c1(:,ikeep),2)); % 1 drop mean spectrum, data 2, condition 1
    Cm12c1(:,i)=abs(Sm12c1./sqrt(Sm1c1.*Sm2c1)); % 1 drop coherence, condition 1
    z1i(:,i)=atanh(Cm12c1(:,i))-1/(dof1-4); % 1 drop, bias-corrected Fisher z, condition 1
    dz1i(:,i)=(z1i(:,i)-z2)/sqrt(1/(dof1-4)+1/(dof2-2)); % 1 drop, z statistic, condition 1
    ps1(:,i)=m1*dz-(m1-1)*dz1i(:,i);
%      ps1(:,i)=dof1*dz-(dof1-2)*dz1i(:,i);
end; 
ps1m=mean(ps1,2);
for j=1:m2;
    jkeep=setdiff(samples2,j); % all samples except j
    Sm12c2=squeeze(mean(S12c2(:,jkeep),2)); % 1 drop mean cross-spectrum, condition 2
    Sm1c2=squeeze(mean(S1c2(:,jkeep),2)); % 1 drop mean spectrum, data 1, condition 2
    Sm2c2=squeeze(mean(S2c2(:,jkeep),2)); % 1 drop mean spectrum, data 2, condition 2
    Cm12c2(:,j)=abs(Sm12c2./sqrt(Sm1c2.*Sm2c2)); % 1 drop coherence, condition 2
    z2j(:,j)=atanh(Cm12c2(:,j))-1/(dof2-4); % 1 drop, bias-corrected Fisher z, condition 2
    dz2j(:,j)=(z1-z2j(:,j))/sqrt(1/(dof1-2)+1/(dof2-4)); % 1 drop, z statistic, condition 2
    ps2(:,j)=m2*dz-(m2-1)*dz2j(:,j);
%      ps2(:,j)=dof2*dz-(dof2-2)*dz2j(:,j);
end;

%
% Leave one out, both samples
% and pseudo values
% for i=1:m1;
%     for j=1:m2;
%         dzij(:,i,j)=(z1i(:,i)-z2j(:,j))/sqrt(1/(dof1-4)+1/(dof2-4));
%         dzpseudoval(:,i,j)=m1*m2*dz-(m1-1)*m2*dz1i(:,i)-m1*(m2-1)*dz2j(:,j)+(m1-1)*(m2-1)*dzij(:,i,j);
% %         dzpseudoval(:,i,j)=dof1*dof2*dz-(dof1-2)*dof2*dz1i(:,i)-dof1*(dof2-2)*dz2j(:,j)+(dof1-2)*(dof2-2)*dzij(:,i,j);
%     end;
% end;
% dzah=sum(sum(dzpseudoval,3),2)/(m1*m2);
ps2m=mean(ps2,2);
% dzar=(sum(ps1,2)+sum(ps2,2))/(m1+m2);
vdz=sum((ps1-ps1m(:,ones(1,m1))).*(ps1-ps1m(:,ones(1,m1))),2)/(m1*(m1-1))+sum((ps2-ps2m(:,ones(1,m2))).*(ps2-ps2m(:,ones(1,m2))),2)/(m2*(m2-1));
% vdzah=sum(sum((dzpseudoval-dzah(:,ones(1,m1),ones(1,m2))).*(dzpseudoval-dzah(:,ones(1,m1),ones(1,m2))),3),2)/(m1*m2);
%
% Test whether H0 is accepted at the specified p value
%
Adz=zeros(size(dz));
x=norminv([p/2 1-p/2],0,1);
indx=find(dz>=x(1) & dz<=x(2)); 
Adz(indx)=1;

if strcmp(plt,'y');
    if isempty(f) || nargin < 6;
        f=linspace(0,1,length(dz));
    end;
    %
    % Compute the coherences
    %
    S121=mean(conj(J1c1).*J2c1,2);
    S122=mean(conj(J1c2).*J2c2,2);
    S111=mean(conj(J1c1).*J1c1,2);
    S221=mean(conj(J2c1).*J2c1,2);
    S112=mean(conj(J1c2).*J1c2,2);
    S222=mean(conj(J2c2).*J2c2,2);
    C121=abs(S121)./sqrt(S111.*S221);
    C122=abs(S122)./sqrt(S112.*S222);
    %
    % Plot the coherence
    %
    subplot(311); 
    plot(f,C121,f,C122); legend('Data 1','Data 2');
    set(gca,'FontName','Times New Roman','Fontsize', 16);
    ylabel('Coherence');
    title('Two group test for coherence');
    subplot(312);
    plot(f,dz);
    set(gca,'FontName','Times New Roman','Fontsize', 16);
    ylabel('Test statistic');
    conf=norminv(1-p/2,0,1);
    line(get(gca,'xlim'),[conf conf]);
    line(get(gca,'xlim'),[-conf -conf]);
    subplot(313);
    plot(f,vdz);
    set(gca,'FontName','Times New Roman','Fontsize', 16);
    xlabel('frequency'); ylabel('Jackknifed variance');
end;
% Adzar=zeros(size(dzar));
% indx=find(dzar>=x(1) & dzar<=x(2)); 
% Adzar(indx)=1;
% 
% Adzah=zeros(size(dzah));
% indx=find(dzah>=x(1) & dzah<=x(2)); 
% Adzah(indx)=1;