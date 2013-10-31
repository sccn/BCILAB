function w = window_func(name,m,param)
% Design a window for a given window function
%
% In:
%   Name : name of the window, can be any of the following:
%          'bartlett' : Bartlett window
%          'barthann' : Bartlett-Hann window
%          'blackman' : Blackman window
%          'blackmanharris' : Blackman-Harris window
%          'flattop'  : Flat-top window
%          'gauss'    : Gaussian window with parameter alpha (default: 2.5)
%          'hamming'  : Hamming window
%          'hann'     : Hann window
%          'kaiser'   : Kaiser window with parameter beta (default: 0.5)
%          'lanczos'  : Lanczos window
%          'nuttall'  : Blackman-Nuttall window
%          'rect'     : Rectangular window
%          'triang'   : Triangular window
%
%   N : number of points in the window
%
%   Param : window parameter (if any)
%
% Out:
%   W : designed window (column vector)
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2013-08-16

% Copyright (C) Christian Kothe, SCCN, 2011, christian@sccn.ucsd.edu
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

p = (0:(m-1))/(m-1);
switch name
    case 'bartlett'
        w = 1 - abs(((0:(m-1)) - (m-1)/2)/((m-1)/2));
    case {'barthann','barthannwin'}
        w = 0.62 - 0.48*abs(p-0.5) - 0.38*cos(2*pi*p);
    case 'blackman'
        w = 0.42-0.5*cos(2*pi*p) + 0.08*cos(4*pi*p);
    case 'blackmanharris'        
        w = 0.35875 - 0.48829*cos(2*pi*p) + 0.14128*cos(4*pi*p) - 0.01168*cos(6*pi*p);
    case {'bohman','bohmanwin'}
        w = (1-abs(p*2-1)).*cos(pi*abs(p*2-1)) + (1/pi)*sin(pi*abs(p*2-1));
    case {'flattop','flattopwin'}
        w = 0.2157 - 0.4163*cos(2*pi*p) + 0.2783*cos(4*pi*p) - 0.0837*cos(6*pi*p) + 0.0060*cos(8*pi*p);
    case {'gauss','gausswin'}
        if nargin < 3
            param = 2.5; end        
        w = exp(-0.5*(param*2*(p-0.5)).^2);
    case 'hamming'
        w = 0.54-0.46*cos(2*pi*p);
    case 'hann'
        w = 0.5-0.5*cos(2*pi*p);
    case 'kaiser'
        if nargin < 3
            param = 0.5; end
        w = besseli(0,param*sqrt(1-(2*p-1).^2))/besseli(0,param);
    case 'lanczos'
        w = sin(pi*(2*p-1))./(pi*(2*p-1)); w(isnan(w)) = 1;
    case {'nuttall','nuttallwin'}
        w = 0.3635819 - 0.4891775*cos(2*pi*p) + 0.1365995*cos(4*pi*p) - 0.0106411*cos(6*pi*p);
    case {'rect','rectwin'}
        w = ones(1,m);
    case 'triang'
        w = 1 - abs(((0:(m-1)) - (m-1)/2)/((m+1)/2));
    otherwise
        % fall back to the Signal Processing toolbox for unknown windows
        if nargin < 3
            w = window(name,m);
        else
            w = window(name,m,param);
        end
end

w = w(:);
