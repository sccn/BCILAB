function H = firlp(N,fp,fs,W)
% Design a low-pass fir filter.
% H = firlp(N,fp,fs,W);
% 
% For W=1, yields results pratically identical to firls when designing a low-pass filter.
%
% In:
%   N  : filter order (must be odd)
%
%   Fp : normalized pass-band frequency edge
%
%   Fs : normalized stop-band frequency edge
%
%   W  : relative weight of the stop-band error
%
% Out:
%   H : the designed FIR filter
%
% Notes:
%   Written after Ivan Selesnick's "Linear-Phase Fir Filter Design By Least Squares" tutorial
%   on http://cnx.org/content/m10577/2.6/
%
% Examples:
%   % design a linear-phase low-pass FIR filter with 10 taps, which transitions within 12-15 Hz
%   % stop-band weight is 1 (as in firls).
%   firlp(10, 12*2/samplerate, 15*2/samplerate, 1);
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

% filter length
M = (N-1)/2;
% calc derivatives
d  = [fp+W*(1-fs) fp*sinc(fp*(1:2*M))-W*fs*sinc(fs*(1:2*M))];
% find least-squares solution
H  = ((toeplitz(d((0:M)+1)) + hankel(d((0:M)+1),d((M:2*M)+1)))/2) \ (fp*sinc(fp*(0:M)'));
% mirror to get a linear-phase FIR representation
H  = [H(M+1:-1:2)/2; H(1); H(2:M+1)/2]';

% sinc function...
function y = sinc(x)
y = ones(size(x));
n = x~=0;
y(n) = sin(pi*x(n))./(pi*x(n));
