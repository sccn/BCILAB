function varargout = exp_eval(x,iters)
% Evaluate the given expression data structure (or short: expression).
% Out-Args = exp_eval(Expression, Iterations)
%
% exp_eval operates on data structures that represent expressions, and evaluates them. The canonical 
% expression data structure is of the form:
%  struct
%    .head: function_handle
%    .parts: cell array
%  
% This structure represents the expression head(parts{:}). Expression data structures can be nested,
% and can have their value attached, if that value is a struct (they can have a "payload", in a
% sense). In this case, they are structured in the form:
%  struct
%    .whatever_fields_belong_to_the_value
%    .tracking: struct
%      .expression: the actual expression data structure
%
% This way, EEG and STUDY data sets can be regular value types and expressions at the same time,
% which allows to annotate EEGLAB data sets with an expression which (uniquely) identifies them, and
% which reproduces them. Expressions are created by some of the expressions/exp_* functions, and by
% functions using the exp_beginfun/exp_endfun contract (i.e., BCILAB filters and data set editing
% components). For this reason, a call like io_quickload('/data/test.set') or
% flt_resample(eeg,[],200) do by default return an unevaluated expression data structure, which must
% be evaluated first (using exp_eval()), to get access to its value. In other words,
% exp_eval(io_quickload('/data/test.set')) gives the EEGLAB data set. This can be viewed as "lazy"
% evaluation. BCILABs framework functions (bci_train, bci_predict, bci_preproc, etc.) do this
% implicitly and automatically.
%
%
% In:
%   Expression : an expression
%   Iterations : maximum number of repeated evaluations applied (default: 1)
%                for real symbolic computations, pass Inf here
%
% Out:
%   Out-Args   : the result of the evaluation, may still be unevaluated if the expression could not
%                be completely evaluated
%
% Notes:
%   For the purposes of exp_eval(), every data structure that can be expressed in MATLAB is an 
%   expression.
%   * Almost all data structures, however, are considered to be 'literal' (and atomic) expressions,
%     such as numbers, strings, function handles, but also recursive data structures, such as cell
%     arrays and structs. These literals cannot be further evaluated. Therefore, almost all MATLAB
%     expressions produce literal expressions, such as, for example, fft(randn(100,1)).
%
%   * A subset of expressions, however, yields by default data that represents unevaluated (or
%     partially evaluated) expressions. Such expressions are formed using either functions that
%     follow the exp_beginfun/exp_endfun protocol, or symbols introduced using exp_symbol or
%     exp_declare. The bulk of these expressions "evaluate" into structures with the two fields
%     'head' and 'parts', where 'head' is a function_handle and 'parts' is a cell array. The purpose
%     of exp_eval() is to fully evaluate these unevaluated expressions, by applying the head to the
%     parts ("calling" the function head with parts as its arguments, in MATLAB terms).
%
%   * To ease interaction with other MATLAB code, some expressions may retain, after having been
%     evaluated (into a MATLAB structure), a subfield that holds the original unevaluated expression
%     that reproduces the data when re-evaluated. These expressions are called "impure" expressions,
%     because they can be viewed both as an unevaluated expression and an evaluated expression. They
%     can thus be processed both by expression-oriented code, and by code that expects values.
%
% Usage Notes:
%   * any MATLAB function can be converted into one that is not immediately evaluated, by adding the
%     two lines "if ~exp_beginfun() return; end", and "exp_endfun;" to its definition, before and
%     after the remaining code.
%   * using exp_symbol, exp_declare_symbols, or @name_of_no_defined_function, arbitrary identifiers
%     can be introduced which are not immediately evaluated (but some of which could be evaluated by
%     exp_eval, given that function code and/or values are defined for them).
%
% Further Reading:
%   A good way to understand symbolic expressions and the operations on them is to read the
%   Mathematica online documentation. The expression system shares some commonalities with
%   Mathematica, with the following considerations:
%    * A subset of functions is implemented with analogous semantics, where SomeFunction is called
%      exp_somefunction in MATLAB
%    * MATLAB has a larger palette of atomic expressions (in addition to String, Number, Symbol, it
%      has Value, which reflects structs, matrices, etc.) In the future, more builtin data structure
%      may be reflected as List, SparseArray, etc.
%    * () are used instead of [], and most operators, such as //., :> cannot be expressed in MATLAB,
%      so that their functional forms must be used.
%    * symbols must be formally introduced before use (with exp_declare_symbols), or expressed
%      either as exp_symbol('x') or as @x
%    * expressions of the form f(a)(b) cannot be interpreted by MATLAB. The workaround is
%      x=f(a);x(b).
%    * expressions are not automatically evaluated; this must be done explicitly with exp_eval().
%      also, the evaluator implemented here is less aggressive, and can be invoked for a limited
%      number of steps.
%    * for the implementation of functions, the exp_beginfun/exp_endfun prologue-epilogue pattern is
%      available.
%
% Examples:
%   % evaluate a dataset-loading expression
%   eeg = exp_eval(io_loadset('/data/projects/test.set'))
%
%   % write up a sequence of filter steps, and finally evalate it; note that, if this particular
%   % sequence had been executed previously, and was cached on disk, the result will be directly 
%   % retrieved from the disk cache
%   eeg = io_loadset('/data/projects/test.set')
%   eeg = flt_resample(eeg,200)
%   eeg = flt_iir(eeg,[0.5 2],'highpass')
%   eeg = exp_eval(eeg);
%
% See also:
%   exp_beginfun, exp_endfun, exp_eval_optimized
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-15

global tracking;
if nargin < 2
    iters = 1; end

varargout = {x};
while iters > 0
    % check if we got a canonically represented expression
    if isfield(x,{'head','parts'})
        chead = char(x.head);
        % we did, check if it's an application of a lambda or builtin function
        if chead(1) == '@' || (exist(char(x.head),'builtin') && ~exist(char(x.head),'file'))
            % if so, need to pre-evaluate the parts
            x.parts = cellfun(@(v)exp_eval(v,iters),x.parts,'UniformOutput',false);
            if ~all(cellfun(@is_evaluated,x.parts))
                % not all of the args could be evaluated: return unevaluated
                return; end
        end
        % apply the function
        if isfield(tracking,'debug') && isfield(tracking.debug,'dbstop_if_error') && tracking.debug.dbstop_if_error && nargout > 0
            % rare case to allow for use of "dbstop if error": call the function without a surrounding try/catch
            [varargout{1:nargout}] = x.head(x.parts{:});
        else
            % use the regular calling mode
            try
                varargout = hlp_wrapresults(x.head,x.parts{:});
            catch e
                if ~(strcmp(e.identifier,'MATLAB:UndefinedFunction') && strcmp(e.stack(1).name,'hlp_wrapresults'))
                    % got a legitimate error: throw it                
                    rethrow(e);
                else
                    % otherwise we just return unevaluated
                    return;
                end
            end
        end
    elseif isa(x,'function_handle')
        str = char(x);
        if str(1) ~= '@'
            % regular function handle
            if is_undefined_function(x)
                % and undefined: do a lookup
                varargout = {hlp_resolve(x,x)}; end
        else
            % lambda function
            try
                % try to obtain the function symbol, if it's a symbolic lambda function
                varargout = {hlp_resolve(get_function_symbol(x),x)};
            catch,end
        end
    end
    
    iters = iters-1;
    
    if iters == 0 || isempty(varargout) || utl_same(x,varargout{1}) || (length(varargout) > 1 && is_evaluated(varargout{1}))
        % evaluation terminates only when the expression no longer changes, or the maximum number of
        % iterations is exceeded (or when multiple evaluated args are being returned and the first
        % one is fully evaluated... - because we want to return these additional outputs)
        break;
    else
        % otherwise we start over...
        x = varargout{1};
    end
end

