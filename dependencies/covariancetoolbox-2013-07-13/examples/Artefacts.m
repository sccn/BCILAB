% Example of artefact detection with the riemannian potato.

% load data
clear all;
load('./data/Artefacts.mat');

Nt = size(data.data,1);
window = 3*data.fs; % 3 second window
overlap = 3*data.fs-data.fs/8;

COV=eeg2cov(data.data',window,overlap);

[C, th] = potato_estimation_iterativ(COV);
[~,d] = potato_detection(COV,C,th);
tcov = (window-overlap)/(data.fs):(window-overlap)/(data.fs):(Nt-window)/data.fs;

t = 1/data.fs:1/data.fs:Nt/data.fs;

figure;
plot(tcov,d);
hold on;
plot(tcov,d>th,'r');
plot(t,data.labels,'k');