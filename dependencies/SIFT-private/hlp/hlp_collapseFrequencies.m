function [collapsed peakidx] = hlp_collapseFrequencies(data,collapsefun,freqidx,tindex,fspacing)
%
% Collapse matrix across frequencies. If tindex is provided, then collapse
% across frequencies for the given timeindex.
%
% Input:
%       data:           col vector or matrix of dimension [nfreqs ntimes].
%                       If matrix, collapse separately for each timepoint
%
%       collapsefun:    The method to use for collapsing across frequencies
%                       none:       no collapse
%                       integrate:  numerical integration using trapz
%                       mean:       average
%                       max:        maximum
%                       absmax:     maximum of absolute value
%                       peak:       1-dimensional peak. Returns 0 if no
%                                   peak is found
%
%       freqidx:        (opt) vector of frequency indices to collapse across
%                       If empty or omitted, use all freqs
%
%       ttindex:        (opt) vector of timepoint(s) or single timepoint
%                       If empty or omitted, use all times
%
%       fspacing:       (opt) Frequency spacing for integration. Default: 1
%
% Out:
%       collapsed:      Row vector or single containing collapsed data for
%                       each time point
%
%       peakidx:        If peak or max is used, return frequency index of
%                       peak location
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%       Theoretical Handbook and User Manual.
%       Available at: http://www.sccn.ucsd.edu/wiki/Sift/
% 
% 
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


% handle defaults
if nargin < 5
    fspacing = 1;   end
if nargin < 4
    tindex  = [];   end
if nargin < 3
    freqidx = [];   end

peakidx = [];

sz    = size(data);
if length(sz)>2
    error('data cannot have more than 2 dimensions');
end

if isempty(freqidx)
    freqidx  = 1:sz(1);    end
if isempty(tindex)
    tindex   = 1:sz(2);    end

% select the desired data range
data = data(freqidx,tindex);
sz   = size(data);

if all(sz==1)
    collapsed = data;
    return;
end

% ensure we have a column vector (we always collapse for each col)
if any(sz==1), data = data(:); end

switch lower(collapsefun)
    case 'none'
        collapsed = data;
    case 'integrate'
        collapsed  = trapz(data)*fspacing;
    case 'mean' 
        collapsed  = mean(data);
    case 'max'
        [collapsed peakidx]  = max(data);
    case 'peak'
         w = warning;
         warning('off','all');
         collapsed   = zeros(1,size(data,2));
         peakidx     = nan(size(collapsed));
         for t=1:size(data,2)
             val = [];
             try, [val idx] = findpeaks(data(:,t),'npeaks',1); catch, end
             
            if isempty(val)
                collapsed(t) = 0;
                peakidx(t)   = nan;
            else
                collapsed(t) = val;
                peakidx(t)   = idx;
            end
         end
        warning(w);
  case 'absmax'
        [collapsed peakidx]  = max(abs(data));
  case 'minmax'
        [collapsed peakidx]  = max(abs(data));
        for k=1:size(peakidx,2)
            collapsed(k) = data(peakidx(k),k); 
        end
end

if isempty(collapsed)
  collapsed = 0;
  peakidx   = nan;
end


