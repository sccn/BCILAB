function [signal,state] = flt_siftpipeline(varargin)
% Implements a basic SIFT pipeline 
% (Preprocessing, Dynamical Modeling, Connectivity estimation)
%
% This filter requires the Source Information Flow Toolbox
% version 1.3b or later (http://sccn.ucsd.edu/wiki/SIFT)
%
% A SIFT object will be returned in signal.CAT
% See hlp_sift_emptyset() for object contents
%
%           Tim Mullen, SCCN/INC/UCSD 2012-14

if ~exp_beginfun('filter'), return; end
declare_properties('name','SIFT', 'experimental',true, 'follows',{'flt_ica','flt_fir','flt_iir','flt_project','flt_sourceLocalize'}, 'cannot_precede',{'flt_clean_settings'},'independent_channels',false, 'independent_trials',false);


signal = arg_extract(varargin,{'signal','Signal','EEG'},[],[]);
signal.srcpot = [];
signal.icaweights = [];
g = arg_define(varargin, ...
        arg_norep({'signal','Signal','EEG'}), ...
        arg({'channels','Channels'}, [], [], 'Cell array of channel names to retain.','type','cellstr','shape','row'), ...
        arg_sub({'preproc','Preprocessing'},{'EEG',signal},@pre_prepData,'Pre-processing options'), ...
        arg_subswitch({'modeling','Modeling'},{'Segmentation VAR' 'EEG',signal}, ...
            hlp_getModelingApproaches, 'Select a modeling approach and define parameters.'), ...
        arg_subtoggle({'selModelOrder','AutoSelectModelOrder'},'off',...
        { ...  % FIXME: for some reason setting this as defaults in modelOrderSelection gens an error: {'plot','off','icselector',{'hq'}}
            arg_sub({'modelOrderSelection','ModelOrderSelection'},{'plot','off','icselector',{'hq'}},@est_selModelOrder,'Model order selection control.','suppress',{'modelingApproach'}), ...
            arg_subswitch({'minimizer','Selector'},{'min'}, ...
            {'elbow', {}, ...
             'min',   {}, ...
             'prctile' ...
             {arg({'prclim','PercentileLimits'},90,[1 100],'Upper percentile limit for order selection. If PercentileLimits = L, This selects the model order, p, for which L% of all sampled windows indicate an optimal model order of p or lower.')} ...
             },'Method for determining optimal model order. If "min", then the optimal model order is the one that minimizes the information criterion. If "elbow" then the optimal order is the elbow of the function of informaton criterion versus model order. If "prctile" then percentile is used.'), ...
        },'Automatically select model order'), ...
        arg_subtoggle({'connectivity','Connectivity'},{'EEG',signal,'MODEL',struct([]),'connmethods',{'S','dDTF08'}},@est_mvarConnectivity,'Select connectivity methods'), ...
        arg_subtoggle({'collapseconn','CollapseConn'},'off',@hlp_collapseConn,'Collapse connectivity across time or freq'), ...
        arg_subtoggle({'validation','Validation'},{'EEG',signal},@est_validateMVAR,'Validate Model fit'), ...
        arg({'printValidation','PrintValidation'},false,[],'Print validation output to console'), ...
        arg({'verb','VerbosityLevel'},0,uint32([0 2]),'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
        arg_nogui({'state','State'}) ...
        );

if g.selModelOrder.arg_selection && length(g.selModelOrder.modelOrderSelection.icselector) > 1
    error('SIFT:onl_siftpipeline:MoreThanOneIC','Only one InformationCriterion can be selected for ModelOrderSelection.'); end
            
state = g.state;
preproc = g.signal;
g      = rmfield(g,'signal');

% Implement basic SIFT pipeline (Preprocessing, Modeling, Connectivity)
% -------------------------------------------------------------------------

% select channels
% -------------------------------------------------------------------------
if ~isempty(g.channels) && ~isequal(g.channels,1:preproc.nbchan) && ~isequal(g.channels,{char(zeros(1,0))}) % TODO: remove last isequal hack
    preproc.chanlocs = preproc.chanlocs(g.channels);
    preproc.data = preproc.data(g.channels,:,:);
    preproc.nbchan = size(preproc.data,1);
end


% buffer data if not epoched
% -------------------------------------------------------------------------
if isempty(preproc.epoch)
    % initialize state if necessary
    if isempty(state)
        for f = utl_timeseries_fields(signal)
            state.(f{1}).buffer = zeros(size(signal.(f{1}),1),0); % buffer of data carried over to next call
            state.(f{1}).next_index = 1;                          % window start index relative to beginning of buffer (1-based)
        end
    end
    % get the desired windowing
    winstep = max(1,round(g.modeling.winstep * preproc.srate));
    winlen = max(1,round(g.modeling.winlen * preproc.srate));
    % get the field that serves as source data (this governs the buffering)
    srcFieldMap = struct('Channels','data','Sources','srcpot','Components','icaact');    
    f = {srcFieldMap.(g.preproc.sigtype.arg_selection)};
    % prepend buffer contents to data (consistently for all time-series fields)
    for tf = utl_timeseries_fields(preproc)
        preproc.(tf{1}) = [state.(tf{1}).buffer preproc.(tf{1})]; end
    % determine the correct window positions for the desired source signal
    g.modeling.winStartIdx = state.(f{1}).next_index:winstep:(size(preproc.(f{1}),2)-winlen+1);
    g.modeling.winStartIdx(g.modeling.winStartIdx<1) = [];
    if ~isempty(g.modeling.winStartIdx)
        % buffer all data that will be needed on the next call
        next_index = max(g.modeling.winStartIdx)+winstep;
        for tf = utl_timeseries_fields(preproc)
            state.(tf{1}).buffer = preproc.(tf{1})(:,next_index:end); end
        % decrement next_index based on how much we at the beginning
        state.(f{1}).next_index = next_index - (size(preproc.(f{1}),2) - size(state.(f{1}).buffer,2));
    else
        % if we cannot make a prediction keep everything buffered
        for tf = utl_timeseries_fields(preproc)
            state.(tf{1}).buffer = preproc.(tf{1}); end
    end
    % consistently update signal meta-data
    preproc.pnts = size(preproc.data,2);
    preproc.xmax = preproc.xmin + (preproc.pnts-1)/preproc.srate;
end

if isempty(preproc.epoch) && isempty(g.modeling.winStartIdx)
    % data is empty, return empty model
    % -------------------------------------------------------------------------
    signal.CAT = hlp_sift_emptyset;
    signal.CAT.Conn = hlp_sift_emptyconn;
    signal.CAT.MODEL = hlp_sift_emptymodel;
else

    % pre-process data
    % -------------------------------------------------------------------------
    preproc = pre_prepData('EEG',preproc,g.preproc,'verb',g.verb,'arg_direct',true);

    % get the m-file name of the function implementing the modeling approach
    % -------------------------------------------------------------------------
    modelingFuncName = hlp_getModelingApproaches('mfileNameOnly',g.modeling.arg_selection);

    % determine optimal model order
    % -------------------------------------------------------------------------
    if g.selModelOrder.arg_selection
        g.selModelOrder.modelOrderSelection.modelingApproach = g.modeling;
        IC = est_selModelOrder('EEG',preproc,g.selModelOrder.modelOrderSelection,'arg_direct',true);
        if iscell(IC), IC = IC{1}; end
        icselector = g.selModelOrder.modelOrderSelection.icselector{1};
        switch g.selModelOrder.minimizer.arg_selection
            case 'min'
                popt = ceil(mean(IC.(icselector).popt));                    
            case 'elbow'
                popt = ceil(mean(IC.(icselector).pelbow));
            case 'prctile'
                popt = ceil(prctile(IC.(icselector).popt,g.selModelOrder.minimizer.prclim));
        end
        g.modeling.morder = popt;
    end

    % fit model
    % -------------------------------------------------------------------------
    preproc.data(~isfinite(preproc.data(:))) = 0;
    preproc.CAT.MODEL = feval(modelingFuncName,'EEG',preproc,g.modeling,'verb',g.verb,'arg_direct',true);

    % perform model validation
    % -------------------------------------------------------------------------
    if g.validation.arg_selection
        [whitestats PCstats stability residualstats] = est_validateMVAR('EEG',preproc,g.validation,'arg_direct',true);
        if isempty(whitestats)
            whitenessCriteria = {};
        else
            whitenessCriteria = g.validation.checkWhiteness.whitenessCriteria;
        end
        if g.printValidation
            vis_validation(whitestats,PCstats,stability,whitenessCriteria);
        end
        % store results
        preproc.CAT.validation.whitestats = whitestats;
        preproc.CAT.validation.PCstats    = PCstats;
        preproc.CAT.validation.stability  = stability;
        preproc.CAT.validation.residualstats = residualstats;
    end

    % calculate connectivity
    % -------------------------------------------------------------------------
    if g.connectivity.arg_selection
        preproc.CAT.Conn = est_mvarConnectivity('EEG',preproc,'MODEL',preproc.CAT.MODEL,g.connectivity,'verb',g.verb,'arg_direct',true);
        if g.collapseconn.arg_selection
            [preproc.CAT.Conn, preproc.CAT.Conn.peak_freqs] = hlp_collapseConn('Conn',preproc.CAT.Conn,g.collapseconn,'arg_direct',true);
        end
    else
        preproc.CAT.Conn = [];
    end

    signal.CAT = preproc.CAT;
end

exp_endfun;