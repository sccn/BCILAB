function [y sigma] = eval_dampedOscillator(f0,tau,fs,variance,ypre)
% evaluate the equation of motion of a damped
% oscillator for a given timepoint.
%
% f0:   fundamental frequency (Hz)
% tau:  damping time (e.g. 100). Larger values --> more sinusoidal
% fs:   samping rate
% x:    previous two samples of the signal (e.g. y(t-1:-1:t-2));
% varnum: the index of the channel/variable being model 
% (e.g. x2 --> varnum = 2)
%
% EXAMPLE: simulate 10000 samples of a stochastically forced 5 Hz damped oscillation
%
% tau = 1000;    % damping time in samples
% f0 = 5;       % fundamental freq in Hz
% srate = 100;  % sampling rate
% 
% y = zeros(1,10000);
% y(1:2) = randn(1,2);
% 
% for n=3:length(y), 
%   y(n) = eval_dampedOscillator(f0,tau,srate,1,y(n-1:-1:n-2)+randn(1,2)); 
% end
% 
% figure; subplot(211); plot(y); subplot(212); pwelch(y);
%
% (C) Tim Mullen, May, 2011. SCCN/INC UCSD

% standard deviation of the signal
sigma = sqrt(1/(1-(exp(-1/tau)^2)));

% equation of motion of damped oscillator
y = 2*exp(-1/tau)*cos(2*pi*f0/fs)*ypre(1) - exp(-2/tau)*ypre(2);

