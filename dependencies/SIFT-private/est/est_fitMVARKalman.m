function [MODEL params] = est_fitMVARKalman(EEG,typeproc,varargin)
%
% Fit multivariate autoregressive model to EEG data using a
% Kalman Filter. This function is a wrapper for the mvaar() function from
% the TSA toolbox [2]. See [1] for additional details on VAR model fitting
% and implementation.
%
% Input:
%
%   EEG                Preprocessed EEG structure.
%   typeproc           Reserved for future use. Use 0
%
% Optional:
%
%     'updatecoeff'        Kalman filter update coefficient {def: 0.001}
%     'updatemode'         Kalman filter noise update mode {def: 1}
%     'winStartIdx'        vector of sample points (start of windows) at which to estimate windowed VAR model
%     'morder'             VAR model order
%     'epochTimeLims'      time range to analyze (sec) where 0 = start of the epoch
%     'verb'               verbosity level (0=no output, 1=text, 2=gui)
%     'timer'              estimate time required to fit model
%       downsampleFactor:   Starting from sample t=max(2,k), store only every k Kalman
%                           coefficients (states, etc) where k=downsampleFactor.
%
% Output:
%
%   MODEL structure with
%       .AR             (numvars x coeffs) matrix of VAR coefficients
%       .PE             (numvars x coeffs) prediction error (noise covariance) coefficients
%       .algorithm      string denoting algorithm used for estimation
%       .modelclass     string denoting model class (here, 'mvar')
%       .Kalman         Kalman state object
%
% See Also: est_fitMVAR(), mvaar()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual. Chapters 3,6.
%     Available at: http://www.sccn.ucsd.edu/wiki/Sift
% [2] A. Schloegl. The time series analysis toolbox for octave and matlab,
%     available online at http://www.dpmi.tu-graz.ac. at/~schloegl/matlab/tsa/
%     and http://cvs.sourceforge.net/ viewcvs.py/octave/octave-forge/extra/tsa/.
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


if nargin<3
    help 'est_fitMVARKalman';
end


var = hlp_mergeVarargin(varargin{:});
g = finputcheck(var, hlp_getDefaultArglist('est'), 'est_fitMVAR','ignore','quiet');
if ischar(g), error(g); end
if isempty(g.epochTimeLims), g.epochTimeLims = [0 EEG.CAT.pnts/EEG.srate]; end
if isempty(g.morder) || length(g.morder)>2, error('invalid entry for field ''morder'''); end
% combine structs, overwriting duplicates of g with
if ~isfield(g,'updatecoeff'), g.updatecoeff = 0.001; end
if ~isfield(g,'updatemode'), g.updatemode = 1; end
if ~isfield(g,'downsampleFactor'), g.downsampleFactor = []; end;
if ~isfield(g,'constraints'), g.constraints = struct([]); end    % constraints.D and constraints.d are constraints of form Dx=d
if ~isfield(g,'storage'), g.storage = {}; end
if ~isfield(g,'Kalman'),  g.Kalman = []; end

if isempty(g.downsampleFactor)
    g.downsampleFactor = round(g.winstep*EEG.srate);
end

g.winstep = g.downsampleFactor/EEG.srate;


%     g = catstruct(g,gvar); clear g2;

if nargout > 1, params = g; end

g.winStartIdx  = max(2,g.downsampleFactor):g.downsampleFactor:size(EEG.CAT.srcdata,2);
numWins   = length(g.winStartIdx);

[AR PE RC]  = deal(cell(1,numWins));

if g.timer
    timeElapsed = nan(1,numWins);
else
    timeElapsed = [];
end
if g.timer, tic; end

% if size(EEG.CAT.srcdata,3)>1
%     if g.verb,
%         fprintf('Kalman filter must operate on single-trial data.\n');
%         fprintf('Concatenating NaN-padded trials to form 2D data matrix.\n');
%     end
%
%     % insert NaN boundary marker between trials and make [nchs x npnts]
%     EEG.CAT.srcdata = nanpad(permute(EEG.CAT.srcdata,[2 1 3]),1);
%     EEG.CAT.srcdata = permute(EEG.CAT.srcdata,[2 1]);
% end

MemoryLengthSamples = -1/log(1-g.updatecoeff);
TimeConstant = MemoryLengthSamples/EEG.srate;
if g.verb,
    fprintf('The effective window length is approximately %0.5f seconds (%d samples)\n',TimeConstant,MemoryLengthSamples);
    fprintf('Your step size is %0.3f ms\n',g.winstep*1000);
    if ~isempty(g.constraints)
        fprintf('Using constraints\n');
    end
end


% Fit MODEL model up to g.morder using Kalman filter

[nchs npnts ntr] = size(EEG.CAT.srcdata);

if isempty(g.storage)
    [VAR MODEL.Q2] = mvaar(permute(EEG.CAT.srcdata,[2 1 3]),g.morder,g.updatecoeff,g.updatemode,g.Kalman,g.verb,g.downsampleFactor,g.constraints);
elseif ismember_bc('residuals',lower(g.storage))
    [VAR MODEL.Q2, MODEL.residuals] = mvaar(permute(EEG.CAT.srcdata,[2 1 3]),g.morder,g.updatecoeff,g.updatemode,g.Kalman,g.verb,g.downsampleFactor,g.constraints);
elseif ismember_bc('kalman',lower(g.storage));
    [VAR MODEL.Q2, MODEL.residuals MODEL.Kalman] = mvaar(permute(EEG.CAT.srcdata,[2 1 3]),g.morder,g.updatecoeff,g.updatemode,g.Kalman,g.verb,g.downsampleFactor,g.constraints);
elseif ismember_bc('state_covariance',lower(g.storage))
    [VAR MODEL.Q2,MODEL.residuals,MODEL.Kalman, MODEL.Q1] = mvaar(permute(EEG.CAT.srcdata,[2 1 3]),g.morder,g.updatecoeff,g.updatemode,g.Kalman,g.verb,g.downsampleFactor,g.constraints);
else
    error('unknown storage option %s',g.storage{1});
end

if g.timer
    time = toc;
end

if g.verb
    h=waitbar(0,'Storing VAR coefficients...');
end

q = 0;
for tr=1:ntr
    % for each trial
    for t=1:size(VAR,1)
        % concatenate VAR coefficients and noise covariance matrices for
        % this trial to end of AR/PE/RC series
        AR{q+t} = reshape(VAR(t,:,tr),nchs*g.morder,nchs).';
        PE{q+t} = MODEL.Q2(:,:,t,tr);
        RC{q+t} = [];
        if g.timer, timeElapsed(t) = time/npnts; end
    end
    q = q + t;
    
    if g.verb
        waitbar(tr/ntr,h);
    end
end

if g.verb, close(h); end

MODEL.AR = AR;
MODEL.PE = PE;
MODEL.RC = RC;

if ntr>1
    % concatenate winStartIdx for each trial into a continuous, increasing sequence
    % e.g. If ntr = 2 and winStartIdx = [1 2 3 4] --> [1 2 3 4 1 2 3 4] --> [1 2 3 4 5 6 7 8]
    g.winStartIdx = bsxfun(@plus,g.winStartIdx,(0:g.winStartIdx(end):ntr*g.winStartIdx(end)-1)')';
    g.winStartIdx = g.winStartIdx(:)';
end

MODEL.winStartTimes = (g.winStartIdx-1)/EEG.srate;
MODEL.winlen = 1/EEG.srate;
MODEL.winstep = g.downsampleFactor/EEG.srate;
MODEL.morder = g.morder;
MODEL.algorithm = 'kalman';
MODEL.modelclass = 'mvar';
MODEL.timeelapsed = timeElapsed;
MODEL.updatecoeff = g.updatecoeff;
MODEL.updatemode = g.updatemode;
MODEL.downsampleFactor = g.downsampleFactor;
MODEL.modelapproach = 'State-Space Modeling';

% MODEL.normalize = g.normalize;

