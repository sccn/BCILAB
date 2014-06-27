function show_fit_quality(x,optmu,optsig,optalpha,optbeta,labels)
hist_bins = 150;
stddev_range = [-5 20];
num_taps = 100;
scale_width = 0.33;

if  length(optmu)<=6
    colors = {'r','g','b','c','m','y'};
else
    map = colormap('jet');
    % trim middle section
    %map(end/3:end*2/3,:) = [];
    for k=1:length(optmu)
        colors{k} = map(1 + floor((k-1)/(length(optmu)-1)*(length(map)-1)),:); end
end


x(x<optmu(1) + stddev_range(1)*optsig(1)) = [];
x(x>optmu(1) + stddev_range(2)*optsig(1)) = [];
hold on;
hist(x,hist_bins);
[hh,bc] = hist(x,hist_bins);
p = [];
for k=1:length(optmu)
    taps = min(x):1/num_taps:max(x);
    if ~isempty(optalpha) && ~isempty(optbeta)
        tmp = exp(-((abs(taps-optmu)./optalpha)).^optbeta)*optbeta/(2*optalpha*gamma(1/optbeta));
    else
        tmp = normpdf(taps,optmu(k),optsig(k));
    end
    scalerange = [optmu(k)-optsig(k)*scale_width,optmu(k)+optsig(k)*scale_width];
    cnt = mean(hh(bc>=scalerange(1) & bc<=scalerange(2)));
    line([optmu(k) optmu(k)],[0 cnt*1.3],'Color',colors{k});
    rescale = cnt/mean(tmp(taps>scalerange(1) & taps < scalerange(2)));
    p(end+1) = plot(taps,tmp*rescale,'Color',colors{k});
end
if exist('labels','var')
    legend(p,labels); end
