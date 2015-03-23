function B = design_fir(N,F,A,nfft,W,odd)
% B = design_fir(N,F,A,nFFT,W)
% Design an FIR filter using the frequency-sampling method.
%
% The frequency response is interpolated cubically between the specified
% frequency points.
%
% In:
%   N : order of the filter
%
%   F : vector of frequencies at which amplitudes shall be defined
%       (starts with 0 and goes up to 1; try to avoid too 
%        sharp transitions)
%
%   A : vector of amplitudes, one value per specified frequency
%
%   nFFT : optionally number of FFT bins to use
%
%   W : optionally the window function to use (default: Hamming)
%
%   Odd : whether to design a filter with odd symmetry (default: false)
%
% Out:
%   B : designed filter kernel
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-08-14

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
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

if nargin < 4 || isempty(nfft)
    nfft = max(8192,2^ceil(log(N)/log(2))); end
if nargin < 5
    W = 0.54 - 0.46*cos(2*pi*(0:N)/N); end
if nargin < 6
    odd = false; end

% calculate interpolated frequency response
F = interp1(round(F*nfft),A,(0:nfft),'pchip');

% set phase & transform into time domain
F = F .* exp(-(0.5*N)*sqrt(-1)*pi*(0:nfft)./nfft);
if odd 
    F = F.*(-i); end %#ok<IJCL>
B = real(ifft([F conj(F(end-1:-1:2))]));

% apply window to kernel
B = B(1:N+1).*W(:)';
