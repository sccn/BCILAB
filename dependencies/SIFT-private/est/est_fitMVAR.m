
function [Model cfg] = est_fitMVAR(varargin)
%
% Fit (adaptive) multivariate autoregressive model to EEG data. See [1] for
% details on VAR model fitting and implementations.
%
%
% Output:
%
%   Model structure with
%       .Model          (numvars x coeffs) matrix of VAR coefficients
%       .PE             (numvars x coeffs) prediction error (noise covariance) coefficients
%       .algorithm      string denoting algorithm used for estimation
%       .modelclass     string denoting model class (here, 'mvar')
%
% See Also: pop_est_fitMVAR(), pop_pre_prepData()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapters 3,6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
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

% pre-cache some frequently used sub-arguments
persistent subargs;
if isempty(subargs)
    subargs.algos = hlp_getMVARalgorithms('defaultNameOnly');
    subargs.all_algos = hlp_getMVARalgorithms;
    subargs.sift_domain = hlp_buildMVARHelpText;
end

verb = arg_extract(varargin,{'verb','VerbosityLevel'},[],2);

haveSigProc = hlp_isToolboxInstalled('Signal Processing Toolbox');
if haveSigProc
    taperFcns = {'rectwin','hamming','hann','bartlett','barthannwin',   ...
                 'blackman','blackmanharris','bohmanwin','chebwin',     ...
                 'flattopwin','gausswin','kaiser','nuttallwin',         ...
                 'parzenwin','taylorwin','tukeywin','triang'};
else
    taperFcns = {'rectwin'};   
end

g = arg_define([0 1],varargin, ...
    arg_norep({'EEG','ALLEEG'},mandatory,[],'EEGLAB dataset'), ...
    arg_subswitch({'algorithm','Algorithm'},subargs.algos,subargs.all_algos,{'Select a model fitting algorithm.',subargs.sift_domain},'cat','Modeling Parameters','suppress',{'ModelOrder','OrderSelector','InitialState'}), ...
    arg({'morder','ModelOrder','modelOrder'},10,[1 Inf],'VAR model order.','cat','Modeling Parameters'), ...
    arg_nogui({'winStartIdx','WindowStartIndices'},[],[],'Starting indices for windows. This is a vector of sample points (start of windows) at which to estimate windowed VAR model','cat','Modeling Parameters'), ...
    arg({'winlen','WindowLength'},0.5,[eps Inf],'Sliding window length (sec)','cat','Modeling Parameters'), ...
    arg({'winstep','WindowStepSize'},0.03,[eps Inf],'Window step size (sec)','cat','Modeling Parameters'), ...
    arg({'taperfcn','TaperFunction'},'rectwin',taperFcns,{'Data tapering (windowing) function', ...
                                                            sprintf(['\n' ...
                                                                     'Each data segment (e.g. all data in a sliding window) will be multipled by the data taper to smooth endpoints towards zero.', ...
                                                                     '\n' ...
                                                                     'Available windows are:\n' ...
                                                                          'rectwin \t - Rectangular window. This is equivalent to no taper.\n', ...
                                                                          'hamming \t - Hamming window.\n', ...
                                                                          'hann \t - Hann window.\n', ...
                                                                          'bartlett \t - Bartlett window.\n', ...
                                                                          'barthannwin \t - Modified Bartlett-Hanning window.\n', ...
                                                                          'blackman \t - Blackman window.\n', ...
                                                                          'blackmanharris- Minimum 4-term Blackman-Harris window.\n', ...
                                                                          'bohmanwin \t - Bohman window.\n', ...
                                                                          'chebwin \t - Chebyshev window.\n', ...
                                                                          'flattopwin \t - Flat Top window.\n', ...
                                                                          'gausswin \t - Gaussian window.\n', ...
                                                                          'kaiser \t - Kaiser window.\n', ...
                                                                          'nuttallwin \t - Nuttall defined minimum 4-term Blackman-Harris window.\n', ...
                                                                          'parzenwin \t - Parzen (de la Valle-Poussin) window.\n', ...
                                                                          'taylorwin \t - Taylor window.\n', ...
                                                                          'tukeywin \t - Tukey window.\n', ...
                                                                          'triang \t - Triangular window.\n' ...
                                                                        ])},'cat','Modeling Parameters'), ...
    arg({'epochTimeLims','EpochTimeLimits'},[],[],'Sub-epoch time limits (sec). This is relative to event time (e.g. [-1 2]). Default is the full epoch time range','cat','Modeling Parameters'), ...
    arg({'prctWinToSample','WindowSamplePercent'},100,[1 100],'Percent of windows to sample','cat','Modeling Parameters'), ...
    arg_subtoggle({'normalize','NormalizeData'},[],@pre_normData,'Z-normalize data within windows. Note this is not recommended for short windows','cat','Window Preprocessing'), ...
    arg_subtoggle({'detrend','Detrend'},{}, ...
    {arg({'method','DetrendingMethod'},'constant',{'linear','constant'},{'Detrend data within each window.', ...
    sprintf(['\n' ...
    'Linear: removes the least-squares fit of a straight line.\n' ...
    'Constant: removes the mean from each trial (centering)' ...
    ]) ...
    } ...
    )},'Detrend or center each time window','cat','Window Preprocessing'), ...
    arg({'timer','Timer'},false,[],'Activate timer. Times are stored in EEG.CAT.Model.timeelapsed'), ...
    arg({'setArgDirectMode','SetArgDirectMode'},true,[],'Set arg_direct mode to true. Can improve speed when number of windows is large. Disable if calling this function repeatedly in a tight loop.'), ...
    arg({'verb','VerbosityLevel'},verb,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );



% commit EEG variable to workspace
[data g] = hlp_splitstruct(g,{'EEG'});
arg_toworkspace(data);
clear data;

% do some error-checking
if isempty(g.epochTimeLims)
    g.epochTimeLims = [EEG.xmin EEG.xmax]; end
if ~(all(g.epochTimeLims>=EEG.xmin) && all(g.epochTimeLims<=EEG.xmax))
    error('Epoch time limits must be within the range [%0.3g %0.3g]',EEG.xmin,EEG.xmax); end
if isempty(g.morder) || length(g.morder)>1
    error('invalid entry for field ''morder.'' Make sure the model order is a single integer.'); end

if nargout > 1, cfg = g; end

%% do some adjustments to parameters

% ensure we are using the right model order
g.algorithm.morder = g.morder;

if rem(g.winstep,1/EEG.srate)
    if g.winstep<1/EEG.srate
        g.winstep = 1/EEG.srate;
    else
        % adjust step size to nearest multiple of sampling interval
        g.winstep = g.winstep-rem(g.winstep,1/EEG.srate);
    end
    
    if g.verb,
        fprintf('Adjusting window step size to nearest multiple of sampling interval\n');
        fprintf('\tstep size is now %0.5g sec\n',g.winstep);
    end
end

if rem(g.winlen,1/EEG.srate)
    if g.winlen<1/EEG.srate
        g.winlen = 1/EEG.srate;
    else
        % adjust window length to nearest multiple of sampling interval
        g.winlen = g.winlen-rem(g.winlen,1/EEG.srate);
    end
    
    if g.verb,
        fprintf('Adjusting window length to nearest multiple of sampling interval\n');
        fprintf('\twindow length is now %0.5g sec\n',g.winlen);
    end
end
if g.winlen > EEG.xmax-EEG.xmin
    g.winlen = EEG.xmax-EEG.xmin;
    if g.verb,
        fprintf('Window length exceeds epoch length. Adjusting window length to match epoch length\n');
        fprintf('\twindow length is now %0.5g sec\n',g.winlen);
    end
end
tidx = getindex(EEG.CAT.times,g.epochTimeLims*1000);
if ~all(isequal(EEG.CAT.times(tidx),g.epochTimeLims*1000))
    
    g.epochTimeLims = EEG.CAT.times(tidx)/1000;
    
    if g.verb
        fprintf('Adjusting epoch time limits to match sampling interval\n');
        fprintf('\tepoch limits are now [%0.5g, %0.5g] sec\n',g.epochTimeLims(1),g.epochTimeLims(2));
    end
end

winLenPnts  = round(g.winlen*EEG.srate); % window size in points
winStepPnts = round(g.winstep*EEG.srate);

if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = tidx(1):winStepPnts:(tidx(2)-winLenPnts)+1;
    %g.winStartIdx  =  round((double(g.epochTimeLims(1):g.winstep:g.epochTimeLims(2)-g.winlen)*EEG.srate)+1;
end

if g.prctWinToSample<100
    % randomly select percentage of windows to work with
    randwin = randperm(length(g.winStartIdx));
    randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
    g.winStartIdx = g.winStartIdx(randwin);
end

numWins   = length(g.winStartIdx);

% construct the data taper
g.taper = window_func(g.taperfcn,winLenPnts).';


%% Prepare data for model fitting

% initialize results arrays
[AR PE RC mu th lambdaOpt]  = deal(cell(1,numWins));

if g.verb==2
    waitbarTitle = sprintf('Fitting VAR[%d] model (%s)...', ...
                        g.morder, ...
                        num2str(g.algorithm.arg_selection));
    
    multiWaitbar(waitbarTitle,'Reset');
    multiWaitbar(waitbarTitle, ...
                 'Color', hlp_getNextUniqueColor, ...
                 'CanCancel','on', ...
                 'CancelFcn',@(a,b) disp('[Cancel requested. Please wait...]'));
end

if g.detrend.arg_selection
    % detrend each window separately
    if g.verb, fprintf('%s detrending each window...\n',firstcaps(g.detrend.method)); end
    for ch=1:EEG.CAT.nbchan
        EEG.CAT.srcdata(ch,:,:) = locdetrend_siftmod(squeeze(EEG.CAT.srcdata(ch,:,:)), ...
            EEG.srate,[g.winlen g.winstep],g.detrend.method);
    end
    if g.verb, fprintf('done.\n'); end
end

if g.timer
    timeElapsed = nan(1,numWins);
else
    timeElapsed = [];
end

%% Main loop: fit MVAR models to each window

algFcnName = hlp_nanocache('algos',10,@hlp_getMVARalgorithms,'mfileNameOnly',g.algorithm.arg_selection);

if g.setArgDirectMode && ~strcmp(algFcnName,'mvar_glADMM') % not necessary for glADMM
    % set the arg_direct flag
    % to improve speed
   g = arg_setdirect(g,true);
end

for t=1:numWins
    
    if g.timer, tic; end
    
    % get data for current window
    data = EEG.CAT.srcdata(:,g.winStartIdx(t):g.winStartIdx(t)+winLenPnts-1,:);
    if g.normalize.arg_selection
        % normalize the data window
        data = pre_normData('data',data,'method',g.normalize.method,'verb',0,'arg_direct',true);
    end

    % execute the model-fitting algorithm
    switch nargout(algFcnName)
        case 2
            [AR{t} PE{t}] = feval(algFcnName, ...
                'data',bsxfun(@times,g.taper,data), ...
                g.algorithm,'arg_direct',true);
        case 3
            [AR{t} PE{t} argsout] = feval(algFcnName, ...
                                   'data',bsxfun(@times,g.taper,data),...
                                    g.algorithm,'arg_direct',true);
            if isstruct(argsout)
                % store contents of argsout fields in cell array at index t
                % e.g. fieldname{t} = argsout.(fieldname)
                fnames = fieldnames(argsout);
                for k=1:length(fnames)
                    eval([fnames{k} '{t}=argsout.(''' fnames{k} ''');']);
                end
            end
        otherwise
            error('SIFT:est_fitMVAR:badAlgArgs','%s must output either 2 or 3 arguments',algFcnName);
    end
        
    if g.verb==2
        drawnow;
        % graphical waitbar
        cancel = multiWaitbar(waitbarTitle,t/numWins);
        if cancel
            if strcmpi('yes',questdlg2( ...
                            'Are you sure you want to cancel?', ...
                            'Model Fitting','Yes','No','No'));
                Model = [];
                multiWaitbar(waitbarTitle,'Close');
                return;
            else
                multiWaitbar(waitbarTitle,'ResetCancel',true);
            end
        end
    end
    
    if g.timer, timeElapsed(t) = toc; end
end

%% Do some cleanup
clear('-regexp','mvar_*')
if g.verb==2
    multiWaitbar(waitbarTitle,'Close'); 
end

%% Construct Model object
Model = hlp_sift_emptymodel;

Model.AR = AR;
Model.PE = PE;
Model.RC = RC;
Model.mu = mu;
Model.th = th;
Model.lambdaOpt = lambdaOpt;
Model.winStartTimes = (g.winStartIdx-1)/EEG.srate;
Model.morder        = g.morder;
Model.winstep       = g.winstep;
Model.winlen        = g.winlen;
Model.algorithm     = g.algorithm.arg_selection;
Model.modelclass    = 'mvar';
Model.timeelapsed   = timeElapsed;
Model.normalize     = g.normalize;
Model.modelapproach = 'Segmentation VAR';
Model.taperFcn      = g.taperfcn;

switch lower(g.algorithm.arg_selection)
    case 'group lasso dal/scsa'
        %     Model.ww = ww;
        Model.lambda = g.algorithm.dal_args.lambda;
    case 'group lasso (admm)'
        Model.lambda = g.algorithm.admm_args.lambda;
        Model.rho    = g.algorithm.admm_args.rho;
        Model.alpha  = g.algorithm.admm_args.alpha;
end
