function [x,window2,window3,f2,f3,t,FS] = radar_signal(x,f,FS)
% Mode 1: 
%   x = radar_signal(N)
%       makes a radar pulse with N time points (for a fixed FS)
%       N is an integer
%
%   [x,win2,win3] = ...
%       returns the windows (with amplitudes builtin) for the 
%       two trapezoidal pulses
%
%   [x,win2,win3,f2,f3,t,FS] = ...
%       also returns the frequencys (in Hz) of the trap. pulses
%       (and the time vector t and sampling rate FS)
%
% Mode 2:
%   x_demod = radar_signal(x_input,f)
%       demodulates x_input at the frequency f


%{
    Realistic Doppler radar:
    10 cm band (i.e. 3 GHz)
    pulse width about 1 us, i.e. 1e-6 s, i.e. 1e3 ns
    pulse repetition rate between 12 kHz and 200 kHz
        i.e. 100 kHz =-> 1e-5 s time interval
%}

if nargin < 2
    N = x;
    % ----------------- Make a signal ---------------------- %
% amplitudes:
A3 = 1;
A2 = -10*A3;
A1 = 1000*A3;

% for plotting:
% A1 = 1; A2 = 10; A3 = 10;

% fix a sampling rate in Hz:
FS = 5e9;       % 5 GHz

T = N/FS; 

% pick an undersampling ratio
% (e.g. for a2i, 8 channels * 50 MHz / 5 GHz
%        = 8/100 )
% M = round( .08 * N );

t = 0:1/FS:(T-1/FS);

% Signal 1: continuous wave
f1 = .8345*1e9;
phase1 = 0;
x = A1*sin( 2*pi*f1*t + phase1 );

% Signal 2: trapezoidal pulse
f2 = 2.33243*1e9;
% window2 = A2*trapWindow(N,.3,.33,.9,.98);

% use units of 1 ns
window2 = pulseTrain( t*1e9, 10,-10, 1e3,40, 1e4, A2 );
x = x + window2 .* sin(2*pi*f2*t);

% Signal 3: small trapezoidal pulse
% f3 = 500*1e6;
% window3 = A3*trapWindow(N,.1,.12,.2,.3);
[window3,multipleWindow3] = pulseTrain( t*1e9, 350, -20, 2e3, 70, 2.2e4, A3 );
nWindows = size(multipleWindow3,1);
% limit frequencies from 200 MHz to 2.4 GHz
f3 = (200 + (2400-200)*rand( nWindows, 1))*1e6;
% I don't want frequencies too close to f1
cnt = 0;
while (any( abs(f3-f1) < 100e6 ) || any( abs(f3-f3) < 50e6 )) && cnt < 100
    cnt = cnt + 1;
    f3 = (200 + (2400-200)*rand( nWindows, 1))*1e6;
end


for n = 1:nWindows
    x = x + multipleWindow3(n,:) .* sin( 2*pi*f3(n)*t );
end

% x = x + window3 .* sin(2*pi*f3*t) ;

else
    % ------------------- Demodulate a signal ----------------- &

    % Make a smooth window -- this helps a lot for edge effects.
    N = length(x);
    nWin = round( N / 10 );
    smoothed = hanning(2*nWin);
    if size(x,1) == 1, smoothed = smoothed.'; end
    x(1:nWin) = smoothed(1:nWin) .* x(1:nWin);
    x(end-nWin+1:end) = smoothed(nWin+1:end) .* x( end-nWin+1:end);
    phase = 1e-14;
    x = myDemod( x, f, 1/FS, 1/50/log10(N), phase );
    
%     x = smooth( x,1000 );
%     x = smooth( x,1000, 'loess' );

end
