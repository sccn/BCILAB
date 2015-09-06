function str = sim_dampedOscillator(f0,tau,fs,varnum)
% return a string corresponding to the equation of motion of a damped
% oscillator. This is meant to be compatible with format used in 
% sim_genVARModelFromEq()
%
% f0:   fundamental frequency (Hz)
% tau:  damping time (e.g. 100). Larger values --> more sinusoidal
% fs:   samping rate
% varnum: the index of the channel/variable being model 
% (e.g. x2 --> varnum = 2)
%
% (C) Tim Mullen, May, 2011. SCCN/INC UCSD

str = sprintf('{2*exp(-1/%f)*cos(2*pi*%f/%f)}*x%d(t-1) + -exp(-2/%f)*x%d(t-2)',tau,f0,fs,varnum,tau,varnum);