function show_fit_quality(x,optmu,optsig)
hist_bins = 50;
stddev_range = [-5 20];
num_taps = 100;
scale_width = 0.33;

x(x<optmu + stddev_range(1)*optsig) = [];
x(x>optmu + stddev_range(2)*optsig) = [];
hold on;
hist(x,hist_bins);
taps = min(x):1/num_taps:max(x);
tmp = normpdf(taps,optmu,optsig);
[hh,bc] = hist(x,hist_bins);
scalerange = [optmu-optsig*scale_width,optmu+optsig*scale_width];
cnt = mean(hh(bc>=scalerange(1) & bc<=scalerange(2)));
line([optmu optmu],[0 cnt*1.3]);
rescale = cnt/mean(tmp(taps>scalerange(1) & taps < scalerange(2)));
plot(taps,tmp*rescale,'r');
