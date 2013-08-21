% Sample usage of the pulseTrain.m function, Stephen Becker

% make a grid of time points
T = 10e-6;   % simulate for a microsecond
N = 2^18;   % number of points (i.e. temporal resolution)
t = linspace(0,T,N);
dT = t(2) - t(1);

% 
ns = 1e-9;
mus = 1e-6;
ms = 1e-3;

offset = 10 *ns;
riseTime = 40 *ns;
duration = 400 *ns;
fallTime = 2*riseTime;
repetitionTime = 1500 *ns;

% ignore the negative in the following (but keep it!)
pulse = pulseTrain( t, offset, -riseTime, duration, fallTime, repetitionTime);

plot(t,pulse,'.');


%{
    Other useful functions from the signal toolbox in Matlab:

type
    help tocwaveformgeneration
to see available functions

%}

%% example, using the chirp (sweeps frequencies)
MHz = 1e6;
GHz = 1e9;

f_initial = 100 * MHz;
f_final = 7.5 * GHz;

% make a grid of time points
T = 10e-0;   % simulate for a second
N = 2^20;   % number of points (i.e. temporal resolution)
t = linspace(0,T,N);
FS = 1/( t(2)-t(1));

myChirp = chirp( t, f_initial, T, f_final );

% plot( t, myChirp, '.' );  % you can't see much
% spectrogram( myChirp, [],[],N/8, FS );