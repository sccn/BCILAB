function x = utl_check_fingerprint(x,opts,ctx,exp)
% Check whether the given argument is an inconsistent impure expression.
% Data = utl_check_fingerprint(Data,Options,Context,Expressions)
%
% The remaining arguments are used only when the function is used as an argstep in exp_beginfun.
% There, it serves as an argstep for 'filter'/'editing' functions, to check, fix up and warn about the value of 
% inconsistent impure expressions (expressions referring to signals, in particular) and data set values 
% (without any notion of expressions).
%
% See also:
%   hlp_fingerprint, exp_beginfun, exp_endfun
%
%                                        Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                        2010-04-15

if ~exist('opts','var')
    opts.fingerprint_check = true;
    exp = [];
    ctx = [];
end

if isfield(x,'tracking') && all(isfield(x.tracking,{'expression','fingerprint'}))
    % get the fingerprint checking expression (either from the options or from the dynamic context) 
    if isempty(opts.fingerprint_check)
        opts.fingerprint_check = hlp_resolve('fingerprint_check',true,ctx); end
    % this check can be selectively overridden with a block-scoped symbol 'fingerprint_check', which may be an expression that involves @expression
    % (the expression for which we are about to execute the check)
    if hlp_microcache('fprint_lookup',@is_enabled,opts.fingerprint_check,exp)
        if ~isequal(hlp_fingerprint(rmfield(x,'tracking')),x.tracking.fingerprint)
            error('expression is inconsistent with the attached value'); end
    end
end

% find out whether a given reference expression yields true if @expression is substituted with some substitution expression
function res = is_enabled(ref_exp,subs_exp)
if exp_eval(utl_releasehold(utl_replacerepeated(ref_exp,{exp_rule(@expression,subs_exp)})),inf)
    res = true;
else
    res = false;
end
