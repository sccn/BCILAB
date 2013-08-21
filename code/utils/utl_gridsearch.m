function [minidx, inputs, outputs] = utl_gridsearch(arg1, varargin)
% Exhaustive search over multiple possible arguments to a function.
% [MinIndex, Inputs, Outputs] = utl_gridsearch(Function, Domain...)
% [MinIndex, Inputs, Outputs] = utl_gridsearch(RangeFormat, Function, Domain...)
% [MinIndex, Inputs, Outputs] = utl_gridsearch(Options, Domain...)
%
% Grid search is a brute-force and assuption-free approach to optimizing parameters of a function.
% It amounts to testing all the possibilities (out of a predetermined set), obtaining the function
% value for each, and returning the parameter combination which gave the best result (here: smallest
% function value). 
%
% The arguments to try for the given function are specified as a regular argument list, in which, 
% depending on the chosen RangeFormat, either arrays among the arguments are interpreted as multiple
% possibilities to try, or arguments which contain a search() expression/data structure are 
% interpreted as defining multiple possibilities. The most convenient way to specify multiple 
% possibilities is by using the special function search(), which packs its multiple arguments into
% an annotated container dat structure which is interpreted by utl_gridsearch as a list of 
% possibilities to try.
%
% Grid search can be run in parallel across multiple processes/machines, by passing the appropriate 
% parallelization options.
%
%
% In:
%   RangeFormat : argument range format, optional; either 'direct' (search ranges are directly 
%                 specified as arrays), or 'clauses' (search ranges are specified using the search()
%                 clause). default: 'direct'
%
%   Function    : the function to optimize; should always return a real number (or NaN) as first 
%                 output, taken as the objective value of the function; may have additional output
%                 arguments, which are also captured and optionally returned by utl_gridsearch
%
%   Options     : cell array of name-value pairs to specify detailed options (alternative to 
%                 RangeFormat/Function); possible names include:
%                 'func'    : the function to optimize (see Function, mandatory) 
%
%                 'argform' : argument range format, either 'direct' or 'clauses' (see RangeFormat, 
%                             default: direct)
%
%                 parallelization options (see par_beginschedule)
%                 'engine_gs' : the parallelization engine to be used for the grid search 
%                               (default: 'local')
%                 'pool'    : node pool to be used for parallelization, when using the BLS scheduler 
%                             (default: 'global')
%                 'policy'  : scheduling policy to be used, when using the BLS scheduler 
%                             (default: 'global')
%
%   Domain...   : one-dimensional argument ranges to be searched; 
%                 * in 'direct' mode, each array determines the possible values at that place in the 
%                   function's argument list, all combinations of specified values are tried; for
%                   cell arrays, the search runs through all given cell contents at that position
%                 * in 'clauses' mode, search ranges for a given argument are specified using the 
%                   search() expression in that place of the argument list
%
% Out:
%   MinIndex    : index of the minimum function return value (first return value if multiple ones 
%                 are returned), out of all tried executions (see example); can be used as index
%                 into Inputs and/or Outputs
%
%   Inputs      : cell array of all tried input combinations (each a cell array of arguments)
%
%   Outputs     : cell array of all received function outputs, one for each input (each a cell array 
%                 of outputs, since the function can have multiple outputs)
%
% Example:
%   % For the four equivalent calls,
%     utl_gridsearch(@f, 1:2, {'a',{'b'},5:7,{5:7}}, [], {}, 'xy', {[1 2 3]});
%     utl_gridsearch('direct', @f, 1:2, {'a',{'b'},5:7,{5:7}}, [], {}, 'xy', {[1 2 3]});
%     utl_gridsearch('clauses', @f, search(1,2), search('a',{'b'},5:7,{5:7}), [], {}, search('x','y'), [1 2 3]);
%     utl_gridsearch({'argform','clauses', 'func',@f}, search(1,2), search('a',{'b'},5:7,{5:7}), [], {}, search('x','y'), [1 2 3]);
%
%   % ... the following 16 = 2x4x1x1x2x1 executions of f are compared:
%       f(1,'a',[],{},'x', [1 2 3])
%       f(1,'a',[],{},'y', [1 2 3])
%       f(1,{'b'},[],{},'x', [1 2 3])
%       f(1,{'b'},[],{},'y', [1 2 3])
%       f(1,5:7,[],{},'x', [1 2 3])
%       f(1,5:7,[],{},'y', [1 2 3])
%       f(1,{5:7},[],{},'x', [1 2 3])
%       f(1,{5:7},[],{},'y', [1 2 3])
%       f(2,'a',[],{},'x', [1 2 3])
%       f(2,'a',[],{},'y', [1 2 3])
%       f(2,{'b'},[],{},'x', [1 2 3])
%       f(2,{'b'},[],{},'y', [1 2 3])
%       f(2,5:7,[],{},'x', [1 2 3])
%       f(2,5:7,[],{},'y', [1 2 3])
%       f(2,{5:7},[],{},'x', [1 2 3])
%       f(2,{5:7},[],{},'y', [1 2 3])
%
% See also:
%   search
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-06

% translate the possible formats for the control arguments into the cell-array-of-options form
if iscell(arg1)
    % utl_gridsearch({Options...}, Domain...)
    args = arg1;
elseif isa(arg1,'function_handle')
    % utl_gridsearch(Function, Domain...)
    args = {'func', arg1};
elseif ischar(arg1)
    % utl_gridsearch(RangeFormat, Function, Domain...)
    args = {'argform',arg1,'func',varargin{1}}; varargin = varargin(2:end);
else
    error('Unexpected control options format.');
end

% read options
opts = hlp_varargin2struct(args, 'func',mandatory, 'argform','direct', 'engine_gs','local', 'pool','global', 'policy','global');

% rewrite arguments if given as clauses
if strcmp(opts.argform,'clauses')
    for i=1:length(varargin)
        varargin{i} = hlp_flattensearch(varargin{i},'cell'); end
elseif ~strcmp(opts.argform,'direct')
    error('Unsupported parameter search mode.');
end

% determine dimensions of the input grid
dims = cellfun('length',varargin);
pitches = [1 cumprod(max(1,dims))];
combos = prod(max(1,dims));

args = varargin;
for i=1:length(args)
    % pre-unpack cell args for singleton arguments
    if iscell(args{i}) && dims(i) == 1
        args{i} = args{i}{1}; end
end
for c=1:combos
    for i=find(dims>1)
        % select current value for non-singleton arguments
        args{i} = varargin{i}(1+floor(mod((c-1)/pitches(i),dims(i))));
        % and unpack cell args
        if iscell(varargin{i})
            args{i} = args{i}{1}; end
    end
    % generate tasks
    tasks{c} = [{@hlp_wrapresults,opts.func},args];
    inputs{c} = args;
end

% schedule tasks
outputs = par_schedule(tasks, 'engine',opts.engine_gs,'pool',opts.pool,'policy',opts.policy);

try
    minidx = argmin(cellfun(@(x)x{1},outputs));
catch
    minidx = NaN;
end
