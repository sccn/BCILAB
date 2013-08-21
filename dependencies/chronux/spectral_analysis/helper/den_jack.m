function [m,ll,ul,llj,ulj]=den_jack(X,family,varargin)
% Function to compute smooth estimates of the mean of x using locfit,
% the corresponding confidence intervals, and jackknife estimates of 
% the confidence intervals
% Usage: [m,ll,ul,llj,ulj]=den_jack(x)
%
% Inputs:
% X: data in the form samples x trials
% family: 'density' or 'reg' for regression
%        If the family is density, the entire input matrix X is considered
%        as data. If the family is regression then the first column of X is
%        taken to be the independent variable and the remaining columns are
%        regressed on this variable (for example, the first column may be
%        the centers of the bins for binned spike count data)
% varargin is the set of arguments used by locfit to perform the smoothing
%
% Outputs:
% m : smoothed estimate of the mean
% ll : estimate of the lower confidence level
% ul : estimate of the upper confidence level
% llj : jackknife estimate of the lower confidence level (+2\sigma
%       where sigma is the jackknife variance)
% llu : jackknife estimate of the upper confidence level (-2\sigma
%       where sigma is the jackknife variance)
[N,NT]=size(X);
if strcmp(family,'reg');
    yy=X(:,2:end);
    y=mean(yy,2);
    x=X(:,1);
    z=scb(x,y,varargin{:});
    figure;
    plot(z(:,1),z(:,2));
    hold on;
    plot(z(:,1),z(:,3),'b:');
    plot(z(:,1),z(:,4),'b:');
    title('Smoothed density estimate, all data');
    
%     fit=locfit(x,y,varargin{:});
%     xfit = lfmarg(fit);
%     yfit = predict(fit,xfit);
%     z = invlink(yfit,fit{4}{5});
% 
    for tr=1:NT-1;
%         i=setdiff(1:NT-1,tr);
%         y=mean(yy(:,i),2);
        y=yy(:,tr);
        fit=locfit(x,y,varargin{:});
        xfit = lfmarg(fit);
        yfit = predict(fit,xfit);
        yfit = invlink(yfit,fit{4}{5});
        zz(:,tr)=yfit;
%         theta(:,tr)=NT*z-(NT-1)*yfit;
    end;    
%     thetam=mean(theta,2);
%     variance=var(theta,0,2);
%     standard_dev=sqrt(variance);
%     figure; plot(xfit{1},thetam,'b');
%     hold on; plot(xfit{1},thetam+2*standard_dev,'r');
%     plot(xfit{1},thetam-2*standard_dev,'r');
%     pause;
    [m,jsd]=jackknife(zz);
%     plot(xfit{1},m,'r');
    hold on; 
    plot(xfit{1},m+2*jsd,'r:');
    plot(xfit{1},m-2*jsd,'r:');
    figure;
    plot(xfit{1},zz);
    title('All trials');
else
    x=mean(X,2);
    fit=locfit(x,varargin{:});
    figure;lfplot(fit);
    lfband(fit);
end;
    


