function W = design_kaiser(lo,hi,atten,odd)
% Design a Kaiser window for a low-pass FIR filter
%
% In:
%   Lo : normalized lower frequency of transition band
%
%   Hi : normalized upper frequency of transition band
%
%   Attenuation : stop-band attenuation in dB (-20log10(ratio))
%
%   OddLength : whether the length shall be odd
%
% Out:
%   W : Designed window
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-08-17

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

% determine beta of the kaiser window
if atten < 21
    beta = 0;
elseif atten <= 50
    beta = 0.5842*(atten-21).^0.4 + 0.07886*(atten-21);
else
    beta = 0.1102*(atten-8.7);
end

% determine the number of points
N = round((atten-7.95)/(2*pi*2.285*(hi-lo)))+1;
if odd && ~mod(N,2)
    N = N+1; end

% design the window
W = window_func('kaiser',N,beta);
