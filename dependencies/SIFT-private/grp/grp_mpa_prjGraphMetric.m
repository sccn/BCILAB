function varargout = grp_mpa_prjGraphMetric(varargin)
% [STUDY obj] = grp_mpa_prjGraphMetric(STUDY,ALLEEG,...)
%
% Perform causal projection [1] using the Measure Projection Toolbox [2]
%
% Output: 
%       STUDY:  modified STUDY with MPT dipoleAndMeasure object in STUDY.measureprojection.sift.object
%       obj:    (optional) dipoleAndMeasure object
% 
% make sure that you have loaded STUDY in which all datasets have SIFT measures
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.9
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
% [2] Bigdely-Shamlo, N, Mullen, T, Kreutz-Delgado,  K, Makeig, S, 
%   Measure Projection Analysis: A Probabilistic  Approach  to EEG Source Comparison and Multi-Subject Inference,
%   NeuroImage (2013), Vol. 72, pp. 287-303, doi:   10.1016/j.neuroimage.2013.01.04
%
% Author: Tim Mullen 10-5-2013, SCCN/INC, UCSD. 
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


% extract some stuff from inputs for arg defaults
ALLEEG = arg_extract(varargin,'ALLEEG',2);

if ~isempty(ALLEEG) && isempty(hlp_checkeegset(ALLEEG(1),{'conn'}));
    Conn            = ALLEEG(1).CAT.Conn;
    ConnNames       = hlp_getConnMethodNames(Conn);
    freqRangeDef    = [Conn.freqs(1) Conn.freqs(end)];
    timeRangeDef    = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
    clear Conn;
else
    ConnNames = {''};
    [freqRangeDef, timeRangeDef] = deal([]);
end

g = arg_define([0 2],varargin, ...
        arg_norep('STUDY',mandatory,[],'STUDY object'),  ...
        arg_norep('ALLEEG',mandatory,[],'ALLEEG array'), ...
        arg({'connmethod','Estimator'},ConnNames{1},ConnNames,'Connectivity estimator to project'), ...
        arg({'timeRange','TimeRange'},timeRangeDef,[],'[Min Max] Time range to project (sec). Leave blank to use all time points','shape','row','type','denserealdouble'), ...
        arg({'freqRange','FrequencyRange'},freqRangeDef,[],'[Min Max] Frequency range to project (Hz). Leave blank to use all frequencies','type','expression','shape','row'), ...
        arg({'collapseTime','CollapseTime'},false,[],'Average across time before projection'), ...
        arg({'collapseFreqs','CollapseFreqs'},false,[],'Integrate across frequencies before projection'), ...
        arg({'centerMeasure','CenterMeasure'},false,[],'Remove mean of each measure. Enabling this may affect interpretability'), ...
        arg({'rejOutsideBrain','RejectDipolesOutsideBrain'},true,[],'Reject dipoles outside the brain'), ...
        arg_sub({'graphMetric','GraphMetric'},{},@hlp_computeGraphMeasure,'Graph metric options','suppress',{'srcNodes','targNodes'}), ...
        arg({'enableCache','EnableCaching'},false,[],'Enable caching of MPT object. This could use a lot of memory, but accelerate subsequent calls to this function'), ...
        arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );
    
% basic error checking
if g.collapseTime && g.collapseFreqs
    error('At most one dimension can be collapsed. Otherwise, there is nothing to smooth!'); end
if isempty(g.STUDY)
    error('Please load a STUDY first'); end
if isempty(g.ALLEEG)
    error('ALLEEG cannot be empty. Did you forget to load a STUDY?'); end
if any(arrayfun(@(x) ~isempty(hlp_checkeegset(x,{'conn'})), g.ALLEEG))
   error(sprintf('You must first run SIFT on all datasets in the STUDY.\nALLEEG.CAT.Conn must be present')); end

waitbarstr = ['Creating MPT Object for ',g.graphMetric.graphMetric];
if g.verb==2
    % reset the color list
    hlp_getNextUniqueColor(hsv(10),'reset');
    multiWaitbar(waitbarstr,'Reset','Color',hlp_getNextUniqueColor);
elseif g.verb==1
    fprintf([waitbarstr '...\n']);
end

% create MPT object
if g.enableCache
    dipoleAndMeasure = hlp_microcache('grp_mpa',@pr.dipoleAndMeasureOfStudySIFT,g.STUDY, g.ALLEEG);
else
    dipoleAndMeasure = pr.dipoleAndMeasureOfStudySIFT(g.STUDY, g.ALLEEG);
end
% reject dipoles outside the brain
if g.rejOutsideBrain
    if g.verb, fprintf('Removing dipoles outside the brain\n'); end
    dipoleAndMeasure = dipoleAndMeasure.createSubsetInRelationToBrain();
end

if g.verb, fprintf('Computing graph measures...\n'); end
for sidx = 1:length(g.ALLEEG)
    % get indices of dipoles for this dataset
    dipidx = dipoleAndMeasure.datasetIdAllConditions==sidx;
    % determine source and target node indices
    g.graphMetric.srcNodes  = dipoleAndMeasure.numberInDataset(dipidx);
    g.graphMetric.targNodes = g.graphMetric.srcNodes;
    % compute graph measure [nodes x freq x time]
    [gm nt nf] = ComputeGraphMeasure(g.ALLEEG(sidx).CAT.Conn, g);
    % subtract mean
    if g.centerMeasure, gm  = bsxfun(@minus, gm, mean(gm,1)); end
    % store vectorized measures
    dipoleAndMeasure.linearizedMeasure(:,dipidx) = gm(:,:)';
    % update progress
    if g.verb==2
        multiWaitbar(waitbarstr,sidx/length(g.ALLEEG));
    elseif g.verb==1, fprintf('.'); end
end

% update some fields
dipoleAndMeasure.time = fastif(length(nt)==1,[],nt*1000);
dipoleAndMeasure.frequency = fastif(length(nf)==1,[],nf);
dipoleAndMeasure.estimator = g.connmethod;
dipoleAndMeasure.measureLabel = g.graphMetric.graphMetric;
dipoleAndMeasure.numberOfMeasureDimensions = ndims(gm)-1;

% store results and prepare outputs
STUDY = g.STUDY;
STUDY.measureprojection.sift.object =  dipoleAndMeasure;
varargout{1} = STUDY;
if nargout>1, varargout{2} = dipoleAndMeasure; end

if g.verb==2
    multiWaitbar(waitbarstr,'Close');
    drawnow;
elseif g.verb==1, fprintf('done\n'); end


% Helper functions
% --------------------------------------------------------------
function [gm newtimes newfreqs] = ComputeGraphMeasure(Conn,g)
% compute graph measures from a Connectivity object

if ~iscell(g.connmethod)
    g.connmethod = {g.connmethod};
end

% collapse connectivity matrices
if g.collapseFreqs
    % collapse across freqs
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethod, ...
        'method',{'freq','net'},'frange',g.freqRange,'freqdim',3,'timedim',4,'verb',0);
else
    % only select subset of freqs
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethod, ...
        'method',{'freq','shrinkonly'}, ...
        'frange',g.freqRange, 'freqdim',3,'timedim',4,'verb',0);
end
if g.collapseTime
    % collapse across time
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethod, ...
        'method',{'time','mean'},'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',0);
else
    % only select subset of times
    Conn = hlp_filterConns(Conn,'connmethods',g.connmethod, ...
        'method',{'time','shrinkonly'}, ...
        'trange',g.timeRange,'freqdim',3,'timedim',4,'verb',0);
end
newtimes = Conn.erWinCenterTimes;
newfreqs = Conn.freqs;
% compute graph measure
gm = hlp_computeGraphMeasure('cmatrix',Conn.(g.connmethod{1}),g.graphMetric);
