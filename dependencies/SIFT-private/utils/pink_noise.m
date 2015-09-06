% function [pn, theta] = phase_noise(num_samp, f0, dbc_per_hz, num_taps)
% 
% This function creates noise with a 1/f spectrum. The noise is then
% phase modulated so that it can be mixed with a signal to simulate
% phase noise in the original signal. The noise is specified in power
% per hertz at frequency f0, relative to the power in the carrier.
% 
% References:
% N. J. Kasdin, "Discrete Simulation of Colored Noise and Stochastic
% Processes and 1/f^a Power Law Noise Generation," _Proceedings of
% the IEEE_, May, 1995.
% Roger L. Freeman, _Reference Manual for Telecommunications
% Engineering_.
% M. Schroeder, _Fractals, Chaos, and Power Laws_.
%
% Input/Output parameters:
% num_samp desired number of output samples
% f0 reference frequency (must be in Hz.)
% dbc_per_hz power per hertz relative to carrier at ref. freq.
% num_taps number of filter taps in AR 1/f filter
% (optional; default = 100)
% 
% pn phase-modulated 1/f process
% theta 1/f process (before phase modulation)
% 

% Jeff Schenck 11/21/95
%
% 1/f noise is produced by passing white noise through a filter. The
% resulting spectrum has the form
% 
% Sx(w) = g^2 / w, (pretend that w is an omega)
% 
% where g is the gain applied to the white noise before filtering. If P
% is the desired power in the 1 Hz. band at w0 = 2pi*f0/fs, and W is the
% 1 Hz. bandwidth in radians, we can write
% 
% P = (1/2pi) Sx(w0) W
% = (1/2pi) g^2/w0 (2pi*1/fs)
% = g^2 / 2pi*f0
% 
% => g = sqrt(2pi*f0*P).
% 
% Notice that the result is *independent* of fs!! Look at it this way:
% if the sampling rate is doubled for a given spectrum, the new w0 is
% half the old w0. For 1/f noise, this means that Sx(w0_new) =
% 2*Sx(w0_old). But a 1 Hz. band is half as large (in radians) than it
% was previously, so the product P is the same, and fs drops out of the
% picture.
% 
% The independence with respect to fs is also an indication of the
% fractal nature of pink noise.
% 
% Note that the phase-modulated noise is itself 1/f if the narrowband
% assumption is valid.

function [pn, theta] = pink_noise(num_samp, f0, dbc_per_hz, num_taps)


% Check input.

if dbc_per_hz >= 0
error('Power per Hz. must be negative.');
elseif f0 <= 0
error('Reference frequency must be positive.');
end

if nargin < 4
num_taps = 100;
end


% Generate white noise. Apply gain for desired dBc/Hz. Warn user
% if gain is too large (gain thresholds have been chosen somewhat
% arbitrarily -- needs work).

gain = sqrt(2*pi * f0 * 10^(dbc_per_hz/10));
wn = gain * randn(1,num_samp);

fprintf('Gain applied to white noise = %f.\n', gain);
if gain >= 1
fprintf('WARNING: Narrowband approximation no longer valid.\n');
elseif gain >= .5
fprintf('WARNING: Narrowband approximation on the verge of collapse.\n');
end


% Generate 1/f AR filter and apply to white noise to produce 1/f
% noise.

a = zeros(1,num_taps);
a(1) = 1;
for ii = 2:num_taps
a(ii) = (ii - 2.5) * a(ii-1) / (ii-1);
end
theta = filter(1,a,wn);


% Phase modulate.

pn = exp(1i*theta);

return;