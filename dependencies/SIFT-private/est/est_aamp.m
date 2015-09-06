
function EEG = est_aamp(varargin)

% Obtain the instantaneous analytic amplitude within a given band of a collection of
% processes in EEG.CAT.srcdata. This uses the hilbert transform and EEGLAB's eegfilt().
%
% See Also: eegfilt(), hilbert(), hlp_aamp()
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
        arg({'ampband','AmplitudePassBand'},[80 150],[],'The [lo hi] pass-band to use for amplitude (Hz)','shape','row'), ...
        arg_subtoggle({'normalize','NormalizeData'},{'method', {'ensemble'}},@pre_normData,'Data normalization. Normalize analytic amplitude trials across time, ensemble, or both','cat','Normalization'), ...
        arg({'plot','PlotData'},false,[],'Plot amplitude envelope'), ...
        arg({'verb','Verbosity'},true,[],'Verbose output') ...
    );

% yeild EEG structure to workspace
data = hlp_splitstruct(g,{'EEG'});
arg_toworkspace(data);
clear data;


%% hilbert transform the data
AA = hlp_aamp(EEG, 'AmplitudePassBand', g.ampband, 'verb', g.verb);

% normalize analytic amplitude, if desired
if g.normalize.arg_selection
    AA = pre_normData('data',AA,g.normalize,'verb',g.verb);
end

% concatenate analytic amplitudes to end of dataset
EEG.CAT.srcdata = cat(1, EEG.CAT.srcdata, AA);

% update EEG datastructure to reflect new AA 'channels'
EEG.CAT.curComponentNames = EEG.CAT.curComponentNames(EEG.CAT.curComps);
EEG.data       = cat(1, EEG.data, AA);
EEG.CAT.curComps = [EEG.CAT.curComps EEG.CAT.curComps+EEG.CAT.nbchan];
EEG.nbchan     = 2*EEG.CAT.nbchan;
EEG.CAT.nbchan = 2*EEG.CAT.nbchan;

N = length(EEG.CAT.curComponentNames);
for k=1:N
    EEG.CAT.curComponentNames{k+N} = [EEG.CAT.curComponentNames{k} '_amp'];
end
% EEG.CAT.curComponentNames = cellstr(num2str(EEG.CAT.curComps'))';
N = length(EEG.chanlocs);
for k=1:N
    EEG.chanlocs(k+N) = EEG.chanlocs(k);
    EEG.chanlocs(k+N).labels = [EEG.chanlocs(k).labels ' amp'];
end

EEG.CAT.datamode = 'AnalyticAmplitude';
EEG.CAT.ampband = g.ampband;

if g.plot
    eegplot(EEG.CAT.srcdata,'srate',EEG.srate,'title',sprintf('AmplitudeModulation [%0.5g %0.5g] Hz',g.ampband(1),g.ampband(2)));
    ax = findobj(gcf,'tag','eegaxis');
    set(ax,'YTickLabel',flipud([EEG.CAT.curComponentNames,' ']'));
end

