function [stats] = est_checkMVARStability(varargin)
%
% Test the stability of a fitted VAR model. See [1-2] for mathematical
% details on testing VAR stability. A stable VAR process is also a
% stationary VAR process [2].
%
% Inputs:
%
%   EEG:        EEGLAB data structure
%   MODEL:      SIFT MODEL structure
%   typeproc:   reserved for future use. Use 0
%
% Optional:
%
%   <Name,Value> pairs containing model fitting parameters. See
%   est_fitMVAR(). Generally, these should be left unspecified.
%
% Outputs:
%
%   stats
%       .stability:  [numwindows x 1] vector of results of stability tests. 1
%                    indicates stable VAR process for that window, 0 indicates
%                    an unstable VAR process.
%
%       .lambda:     [numwindows x nchs*morder] matrix of eigenvalues of VAR
%                    process. All eigenvalues should be < 1 for stable VAR process
%
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual. Chapters 3,6.
%     Available at: http://www.sccn.ucsd.edu/wiki/Sift
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
%     Springer.
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD
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

g = arg_define([0 2],varargin, ...
        arg_norep({'EEG','ALLEEG'},mandatory,[],'EEGLAB dataset'), ...
        arg_norep({'MODEL','Model'},mandatory,[],'MVAR MODEL object'), ...
        arg_nogui({'winStartIdx','WindowStartIndices'},[],[],'Starting indices for windows. This is a vector of sample points (start of windows) at which to estimate windowed VAR model. Default is empty (use all windows)','cat','Options'), ...
        arg({'prctWinToSample','WindowSamplePercent'},100,[1 100],'Percent of windows to sample','cat','Options'), ...
        arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
        );
    
% commit EEG and MODEL variables to workspace
[data g] = hlp_splitstruct(g,{'EEG','MODEL'});
arg_toworkspace(data);
clear data;

morder = MODEL.morder;

% window size in points
% winLenPnts = floor(MODEL.winlen*EEG.srate);

if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = round(MODEL.winStartTimes*EEG.srate)+1;
end

if g.prctWinToSample<100 
    % randomly select percentage of windows to work with
    randwin = randperm(length(g.winStartIdx));
    randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
    g.winStartIdx = g.winStartIdx(randwin);
    g.winArrayIndex = randwin;
end

% get the array indices of the windows we are working with
g.winArrayIndex = getindex(MODEL.winStartTimes,(g.winStartIdx-1)/EEG.srate);

% initialize waitbar
if g.verb==2
    waitbarTitle = sprintf('Checking stability %s...', ...
        fastif(isempty(EEG.condition),'',['for ' EEG.condition]));
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle,'ResetCancel',true);
    multiWaitbar(waitbarTitle, ...
                 'Color', hlp_getNextUniqueColor, ...
                 'CanCancel','on', ...
                 'CancelFcn',@(a,b) disp('[Cancel requested. Please wait...]'));
end

numWins = length(g.winStartIdx);

stats.stability = zeros(1,numWins);
[nchs Mp] = size(MODEL.AR{1});
stats.lambda = zeros(numWins,Mp);
%lambda = [];
I = eye(nchs*morder-nchs,nchs*morder-nchs);
O = zeros(nchs*morder-nchs,nchs);
for t=1:numWins
    % get the array index of the window we are working with
    winArrIdx = g.winArrayIndex(t);
        
    % rewrite VAR[p] process as VAR[1]
    A = [MODEL.AR{winArrIdx} ; [I O]];
    stats.lambda(t,:) = log(abs(eig(A)));
    stats.stability(t) = all(stats.lambda(t,:)<0);
    
    if g.verb==2
        % update waitbar
        drawnow;
        cancel = multiWaitbar(waitbarTitle,t/numWins);
        if cancel && hlp_confirmWaitbarCancel(waitbarTitle)
            stats = [];
            return;
        end
    end
    
end

stats.winStartIdx = g.winStartIdx;
stats.winStartTimes = MODEL.winStartTimes(g.winArrayIndex);
stats.winArrayIndex = g.winArrayIndex;

% clean up
if g.verb==2
    multiWaitbar(waitbarTitle,'Close'); 
end


