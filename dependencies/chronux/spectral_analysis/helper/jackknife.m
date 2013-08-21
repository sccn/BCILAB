function [m,jsd]=jackknife(x)
% Compute jackknife estimates of the mean and standard deviation of input data x
% Usage: [m,jsd]=jackknife(x)
% 
% Inputs:
% x : data in the form samples x trials
%
% Outputs:
% m : estimate of the mean (across trials)
% jsd: jackknife estimate of the standard deviation (across trials)

[N,C]=size(x);
if C==1; error('Need multiple trials'); end;
m=mean(x,2);
theta=zeros(N,C);
for tr=1:C;
    i=setdiff((1:C),tr); % drop 1 trial
    y=sum(x(:,i),2)/(C-1); % mean over remaining trials
    theta(:,tr)=C*m-(C-1)*y; % pseudo values
%     yy(:,tr)=y;
end;
jm=mean(theta,2);
jm=repmat(jm,[1 C]);
% jm2=mean(yy,2);
% jm2=repmat(jm2,[1 C]);
jsd=sqrt(sum((theta-jm).^2,2)/(C*(C-1)));
% jsd2=sqrt((C-1)*sum((yy-jm2).^2,2)/C);
% jsd
% jsd2
