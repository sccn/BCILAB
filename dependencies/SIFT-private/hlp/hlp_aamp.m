function [AA] = hlp_aamp(varargin)

% Obtain the instantaneous analytic amplitude within a given freq band for a collection of
% processes in EEG.CAT.srcdata. This uses the hilbert transform and EEGLAB's eegfilt().
%
% Input:
%
%   EEG                 Preprocessed EEG structure.
%   AmplitudePassBand   The [lo hi] pass-band for which to estimate instantaneous amplitude
%
% Outputs
%   AA                  analytic amplitude of same dimensions as
%                       EEG.CAT.srcdata
%
% See Also: eegfilt(), hilbert(), est_PMGC()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/SIFT/
% 
% 
% Author: Tim Mullen 2012, SCCN/INC, UCSD. 
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

g = arg_define([0 1], varargin, ...
        arg_norep({'EEG'},mandatory),...
        arg({'ampband','AmplitudePassBand'},[10 30],[],'The [lo hi] pass-band to use for amplitude (Hz)','shape','row'), ...
        arg({'verb','Verbosity'},2,{int32(0) int32(1) int32(2)},'Verbosity level.') ...
    );


[nchs npnts ntr] = size(g.EEG.CAT.srcdata);

if g.ampband(2)>g.EEG.srate/2
    error('SIFT:est_aamp','Upper pass-band cutoff must be less than the Nyquist rate of %0.5g Hz',g.EEG.srate/2);
end

if g.verb
    fprintf('Computing [%0.2g-%0.2g]Hz analytic amplitude using Hilbert transform...\n',g.ampband(1),g.ampband(2));
end

if g.verb==2
    multiWaitbar('Computing analytic amplitude','Reset','Color',hlp_getNextUniqueColor);
end

AA = zeros(nchs, npnts, ntr);
filtd = eegfilt(g.EEG.CAT.srcdata(:,:),g.EEG.srate,g.ampband(1),[],npnts);
filtd = eegfilt(filtd,g.EEG.srate,[],g.ampband(2),npnts);
filtd = reshape(filtd,[nchs, npnts, ntr]);

if nchs>ntr
    % hilbert transform each trial separately (all channels at once)
    for tr=1:ntr
        if g.verb==2
            multiWaitbar('Computing analytic amplitude',tr/ntr);
        end
        AA(:,:,tr)=abs(hilbert(filtd(:,:,tr)'))'; 
    end
else
    % hilbert transform each channel separately (all trials at once)
    for ch=1:nchs
        if g.verb==2
            multiWaitbar('Computing analytic amplitude',ch/nchs);
        end
        AA(ch,:,:) = abs(hilbert(filtd(ch,:,:))); 
    end
end

if g.verb==2
    multiWaitbar('Computing analytic amplitude','Close');
end
