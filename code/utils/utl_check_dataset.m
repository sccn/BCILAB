function sig = utl_check_dataset(sig,opts,ctx,exp)
% Check whether the given argument is an imporperly tracked data set and fix.
% Data = utl_check_dataset(Data,Options,Context,Expressions)
%
% The remaining arguments are used only when the function is used as an argstep in exp_beginfun.
% There, it serves as an argstep for 'filter'/'editing' functions, to check, fix up and warn about
% the value of inconsistent impure expressions (expressions referring to signals, in particular) and
% data set values (without any notion of expressions).
%
% See also:
%   exp_beginfun, exp_endfun
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-15

if nargin < 2
    opts.fingerprint_check = true;
    ctx = [];
    exp = [];
end

% check if we have a data set
if isfield(sig,{'data','srate'})
    if ~isfield(sig,'tracking') || ~isfield(sig.tracking,'expression')
        % the data set was imported from EEGLAB or the tracking info was ruined
    elseif isfield(sig.tracking,'fingerprint')
        % get the fingerprint checking expression (either from the options or from the dynamic scope)
        if isempty(opts.fingerprint_check)
            opts.fingerprint_check = hlp_resolve('fingerprint_check',true,ctx); end
        % check whether it is enabled
        if ~opts.fingerprint_check
            return; end
        if isequal(opts.fingerprint_check,true) || hlp_microcache('fprint_lookup',@is_enabled,opts.fingerprint_check,exp)
            fprint = hlp_fingerprint(rmfield(sig,'tracking'));
            if fprint ~= sig.tracking.fingerprint
                % the data set has been edited according to some fingerprint check
            else
                % no problem was spotted
                return; 
            end
        else
            % fingerprinting disabled (nothing to check here...)
            return; 
        end
    else
        return; % the data set does hot have a fingerprint
    end
    
    % --- the data set has been imported from EEGLAB or has been edited with it: do an implict flush-to-disk & re-import ---
    
    if isfield(sig,'tracking')
        if isfield(sig.tracking,'expression')
            sig.tracking = rmfield(sig.tracking,'expression'); end
        if isfield(sig.tracking,'online_expression')
            sig.tracking = rmfield(sig.tracking,'online_expression'); end
        if isfield(sig.tracking,'fingerprint')
            sig.tracking = rmfield(sig.tracking,'fingerprint'); end
    else
        sig.tracking = struct();
    end
    
    % check for problems by scanning the history
    found_problem = 0;
    if isfield(sig,'history')
        % catch a few known to be problematic operations in the imported data set
        operations = { ...
            {{'pop_epoch'}, ...
            'Note: The BCILAB framework cannot operate on datasets that have been epoched using pop_epoch in EEGLAB.\n'}, ...
            {{'pop_eegfilt','pop_runica'}, ...
            'Note: The operation ''%s'' will neither allow to generate reliable predictions nor usable online models from the imported data set.\n'}, ...
            {{'pop_averef','pop_interp','pop_reref','pop_resample','pop_rmbase','pop_subcomp'}, ...
            'Note: The operation ''%s'' will not allow to generate online models from the imported data set.\n'}, ...
            {{'pop_chansel','pop_rejchan','pop_rejchanspec','pop_select'}, ...
            ['Note: The operation ''%s'' requires that the channels supplied during online processing (and their order) match those' ...
            ' in the present data set.\n']}, ...
            };
        for l=1:length(operations)
            for op=operations{l}{1}
                if strfind(sig.history,[op{1} '('])
                    fprintf(operations{l}{2},op{1});
                    found_problem = 1;
                end
            end
        end
        if found_problem
            fprintf(['\nIt is recommended that these operations be executed from BCILAB''s ' ...
                'processing palette, since the operations implemented there are online-capable.\n']);
        end
    end
    
    % (re-) create the fingerprint, if necessary
    if ~exist('fprint','var')
        fprint = hlp_fingerprint(rmfield(sig,'tracking')); end
    
    % flush data to disk, if not already there...
    filepath = ['temp:/flushedsets/' num2str(fprint) '.mat'];
    if ~exist(env_translatepath(filepath),'file')
        disp('Flushing data set to disk...');
        EEG = sig; %#ok<NASGU>
        io_save(filepath,'EEG','-makedirs','-attributes','''+w'',''a'''); 
    end
    % change the expression into a re-loading
    sig = io_loadset(filepath);
elseif isfield(sig,{'head','parts'})
    % we have an expression: check parts recursively...
    for p=1:length(sig.parts)
        sig.parts{p} = utl_check_dataset(sig.parts{p},opts,ctx,exp); end
end


% find out whether a given reference expression yields true if @expression is substituted with some substitution expression
function res = is_enabled(ref_exp,subs_exp)
if exp_eval(utl_releasehold(utl_replacerepeated(ref_exp,{exp_rule(@expression,subs_exp)})),inf)
    res = true;
else
    res = false;
end
