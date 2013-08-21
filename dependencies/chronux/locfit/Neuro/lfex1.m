% Local Estimation a spike firing rate (in spikes per unit time).
%
% The estimation procedure approximates the log of the spike firing
% rate by a quadratic polynomial, within a sliding window.
%
% The first argument to locfit() is a column vector of the spike times.
% 'family','rate' specifies that the output of the density estimate should
% be in terms of events per unit time. (instead of the default, statistical
% density estimation).
%
% 'alpha',0.6 specifies the width of the sliding windows, as a fraction of
% the total spikes. 
%
% xlim gives the limits of the observation interval. Correct specification
% of this is critical to avoid introducing bias at end-points.
%
% The lfband(fit) line adds confidence bands (based on pointwise 95%
% coverage) to the plot.
%
% The spike time data is from Hemant Bokil.

load lmem5312.mat;
fit = locfit(data(6).times{1},'xlim',[78.9 80.30],'family','rate','nn',0.6);
lfplot(fit);
title('Spike Firing Rate Estimation');
lfband(fit);