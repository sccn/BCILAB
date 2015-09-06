function [Conn g] = est_mvarConnectivity(varargin)
%
% Calculate spectral, coherence, and connectivity measures from a fitted
% VAR model. See [1] for additional details on VAR model fitting and
% connectivity measures implemented here.
%
% Inputs:
%
%   EEG:        EEGLAB data structure
%   MODEL:      VAR model returned by est_fitMVAR or other SIFT model-fitting
%               routine
%   connmethods: A cell array of connectivity method names compatible with
%               the names in est_mvtransfer().
%               
%
% Optional:
%
%   <Name,Value> pairs containing model fitting parameters. See
%   est_fitMVAR(). Generally, these should be left unspecified and inferred
%   from MODEL structure.
%
% Output:
%
%   Conn
%       .<measure>        [numvar x numvar x numfreq x numtime] connectivity
%                         matrix
%       .freqs            vector of frequencies estimated
%       .erWinCenterTimes vector of window centers (seconds) relative to
%                         event onset
%   params                parameters used
%
% See Also: pop_est_mvarConnectivity(), pop_est_fitMVAR(), est_mvtransfer()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%     Theoretical Handbook and User Manual. Chapter 4
%     Available at: http://www.sccn.ucsd.edu/wiki/Sift
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

if ~isempty(varargin)
    EEG = arg_extract(varargin,{'EEG', 'ALLEEG'},1);
else
    EEG = [];
end

if ~isempty(EEG) && isfield(EEG,'srate')
    if EEG.srate > 1
        defreqs = ['1:1:' num2str(fix(EEG.srate/2)-1)];
    else
        defreqs = ['0.01:0.01:' num2str(EEG.srate/2-0.01)];
    end
else
    defreqs = [];
end

% arg_sub({'connmethods','ConnectivityMeasures','ConnectivityMethods'},[],@est_mvtransfer,'Select measures.','suppress',{'freqs'},'cat','Options'), ...

g = arg_define([0 2],varargin, ...
        arg_norep({'EEG','ALLEEG'},mandatory,[],'EEGLAB dataset','cat','Data Input'), ...
        arg_norep({'MODEL','Model'},mandatory,[],'MVAR MODEL object','cat','Data Input'), ...
        arg({'connmethods','ConnectivityMeasures','ConnectivityMethods'},{},hlp_getValidConnMethods, ...
                    {'Select measures to estimate', ...
                    sprintf(['\n' ...
                             '\n' ...
                             'Measures are categorized as follows:\n' ...
                                '+ DIRECTED TRANSFER FUNCTION MEASURES\n', ...
                                '     DTF:    Directed Tranfer Function\n',...
                                '     nDTF:   Normalized DTF\n',...
                                '     dDTF:   Direct DTF\n',...
                                '     dDTF08: Direct DTF (with full causal normalization)\n',...
                                '     ffDTF:  Full-frequency DTF\n',...
                                '+ PARTIAL DIRECTED COHERENCE MEASURES\n', ...
                                '     PDC:    Partial Directed Coherence\n',...
                                '     nPDC:   Normalized PDC\n', ...
                                '     GPDC:   Generalized Partial Directed Coherence\n', ...
                                '     PDCF:   Partial Directed Coherence Factor\n',...
                                '     RPDC:   Renormalized Partial Directed Coherence\n', ...
                                '+ GRANGER-GEWEKE CAUSALITY MEASURES\n', ...
                                '     GGC:    Granger-Geweke Causality\n', ...
                                '+ SPECTRAL COHERENCE MEASURES\n',...
                                '     Coh:    Complex Coherence\n', ...
                                '     iCoh:   Imaginary Coherence\n', ...
                                '     pCoh:   Partial Coherence\n', ...
                                '     mCoh:   Multiple Coherence\n', ...
                                '+ SPECTRAL DENSITY MEASURES\n', ...
                                '     S:      Complex Spectral Density\n' ...
                                ])}, ...
                        'type','logical','cat','Connectivity Estimation'), ...
        arg({'absvalsq','SquaredModulus','AbsValSquared'},true,[],'Square the modulus. Return the squared modulus (square of the absolute value) of complex measures.','cat','Options'), ...
        arg({'spectraldecibels','ConvertSpectrumToDecibels'},false,[],'Return spectral power in decibels','cat','Options'), ...
        arg({'freqs','Frequencies'},defreqs,[],'Frequencies. All measures will be estimated and returned for each of these frequencies','type','expression','shape','row','cat','Options'), ...
        arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
        );
    
% commit EEG and MODEL variables to workspace
[data g] = hlp_splitstruct(g,{'EEG','MODEL'});
arg_toworkspace(data);
clear data;

% do some input checking
if isempty(MODEL)
    error('SIFT:est_mvarConnectivity','You must fit an MVAR model before estimating connectivity. See pop_est_fitMVAR().');
end
if isempty(g.connmethods)
    error('SIFT:est_mvarConnectivity','Please select at least one connectivity measure');
end
if any(g.freqs<0) || any(g.freqs>EEG.srate/2)
    error('SIFT:est_mvarConnectivity','Frequencies must be within the range [%d %d] Hz',0,EEG.srate/2);
end

g.winstep   = MODEL.winstep;
g.winlen    = MODEL.winlen;

if isempty(g.freqs), g.freqs = eval(defreqs); end

numWins = length(MODEL.AR);
nchs    = EEG.CAT.nbchan;

% initialize Conn object
Conn = hlp_sift_emptyconn;

if g.verb
    % inform user of total amount of memory (megabytes) required
    bytesReq = 4*length(g.connmethods)*numWins*length(g.freqs)*nchs^2;
    
    fprintf('-------------------------------------------------------------------------\n');
    fprintf(['Connectivity estimation will require %5.5g MB of memory (per condition).\n' ...
             'Make sure you have enough memory available.\n'],bytesReq/(1024^2));
    fprintf('-------------------------------------------------------------------------\n');

    % check if we are likely to exceed memory capacity and notify user
    bytesAvail = hlp_getAvailableMemory('bytes');
    if bytesReq > bytesAvail
        if g.verb==2
            res = questdlg2('It appears you may not have sufficient memory to carry out this operation. Do you want to continue?','est_mvarConnectivity: Memory check','Yes','No','No');
        else
            res = input('It appears you may not have sufficient memory to carry out this operation.\nDo you want to continue? (y)es/(n)o:  ','s');
        end
        
        if strcmpi(res(1),'n')
            return;
        end
    end
end

if g.verb==2
    % create waitbar
    waitbarTitle = sprintf('Estimating connectivity %s...', ...
        fastif(isempty(EEG.condition),'',['for ' EEG.condition]));
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle, ...
                 'Color', hlp_getNextUniqueColor, ...
                 'CanCancel','on', ...
                 'CancelFcn',@(a,b) disp('[Cancel requested. Please wait...]'));
end
    
    
for t=1:numWins
    
    if isempty(MODEL.AR{t})
        continue;
    end
    
    % extract noise covariance matrix
    if size(MODEL.PE{t},2)>nchs
        MODEL.PE{t} = MODEL.PE{t}(:,nchs*MODEL.morder+1:nchs*(MODEL.morder+1));
    end
    
    % estimate connectivity for this window
    ConnTmp = est_mvtransfer('AR',MODEL.AR{t}, ...
                             'C',MODEL.PE{t},  ...
                             'freqs',g.freqs,  ...
                             'srate',EEG.srate,...
                             'connmethods', g.connmethods, ...
                             'arg_direct', true);
    
    for method=g.connmethods
        if t==1
            % on first run, preallocate memory for connectivity matrices
            try
                Conn.(method{1}) = zeros([size(ConnTmp.(method{1})) numWins],'single');
            catch
                errordlg2('Insuficient memory available to allocate variables');
                return;
            end
        end
        
        Conn.(method{1})(:,:,:,t) = ConnTmp.(method{1});
    end
    
    
    if g.verb==2
        % graphical waitbar
        drawnow;
        cancel = multiWaitbar(waitbarTitle,t/numWins);
        if cancel && hlp_confirmWaitbarCancel(waitbarTitle)
            Conn = [];
            return;
        end
    end
    
end

% transform the connectivity as user requested
if g.absvalsq
    if g.verb
        fprintf('Returning squared modulus of complex measures.\n'); 
    end
    Conn = hlp_absvalsq(Conn,g.connmethods,false,g.verb);
    % special case for iCoh
    Conn = hlp_absvalsq(Conn,{'iCoh'},true,0);
end

if g.spectraldecibels && isfield(Conn,'S')
    if g.verb
        fprintf('Returning spectrum in decibels.\n'); 
    end
    Conn.S = 10*log10(Conn.S);
end

if g.verb
    fprintf('Creating final connectivity object %s...\n', ...
            fastif(isempty(EEG.condition),'',['for ' EEG.condition]));
end

if ~strcmpi(MODEL.algorithm,'kalman')
    Conn.winCenterTimes = MODEL.winStartTimes+MODEL.winlen/2;
else
    Conn.winCenterTimes = MODEL.winStartTimes;
end

Conn.erWinCenterTimes = Conn.winCenterTimes-abs(EEG.xmin);
Conn.freqs = g.freqs;
    
if g.verb==2
    multiWaitbar(waitbarTitle,'Close'); 
end

