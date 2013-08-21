function res = hlp_aggregatestructs(structs,defaultop,varargin)
% Aggregate structs (recursively), using the given combiner operations.
% Result = hlp_aggregatestructs(Structs,Default-Op,Field-Ops...)
%
% This results in a single 1x1 struct which has aggregated values in its fields (e.g., arrays, 
% averages, etc.). For a different use case, see hlp_superimposedata.
%
% In:
%   Structs    : cell array of structs to be aggregated (recursively) into a single struct
%
%   Default-Op : optional default combiner operation to execute for every field that is not itself a 
%                struct; see notes for the format.
%
%   Field-Ops  : name-value pairs of field-specific ops; names can have dots to denote operations
%                that apply to subfields. field-specific ops that apply to fields that are
%                themselves structures become the default op for that sub-structure
%
% Out:
%   recursively merged structure.
%
% Notes:
%   If an operation cannot be applied, a sequence of fall-backs is silently applied. First,
%   concatenation is tried, then, replacement is tried (which never fails). Therefore,
%   function_handles are being concatenated up to 2008a, and replaced starting with 2008b.
%   Operations are specified in one of the following formats:
%   * 'cat': concatenate values horizontally using []
%   * 'replace': replace values by those of later structs (noncommutative)
%   * 'sum': sum up values
%   * 'mean': compute the mean value
%   * 'std': compute the standard deviation
%   * 'median': compute the median value
%   * 'random': pick a random value
%   * 'fillblanks': replace [] by values of later structs
%   * binary function: apply the function to aggregate pairs of values; applied in this order
%     f(f(f(first,second),third),fourth)...
%   * cell array of binary and unary function: apply the binary function to aggregate pairs of
%     values, then apply the unary function to finalize the result: functions {b,u} are applied in
%     the following order: u(b(b(b(first,second),third),fourth))
%
% Examples:
%   % calc the average of the respective field values, across structs
%   hlp_aggregatestructs({result1,result2,result3},'mean')
%
%   % calc the std deviation of the respective field values, across structs
%   hlp_aggregatestructs({result1,result2,result3},'std')
%
%   % concatenate the field values across structs
%   hlp_aggregatestructs({result1,result2,result3},'cat')
%
%   % as before, but use different operations for a few fields
%   hlp_aggregatestructs({result1,result2,result3},'cat','myfield1','mean','myfield2.subfield','median')
%   
%   % use a custom combiner operation (here: product)
%   hlp_aggregatestructs({result1,result2,result3},@times)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-04

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

warning off MATLAB:warn_r14_function_handle_transition

if ~exist('defaultop','var')
    defaultop = 'cat'; end
if ~iscell(structs)
    structs = {structs}; end
fieldops = hlp_varargin2struct(varargin);

% translate all ops (if they are specified as strings)
defaultop = translateop(defaultop);
fieldops = translateop(fieldops);

% aggregate, then finalize
res = finalize(aggregate(structs,defaultop,fieldops),defaultop,fieldops);


function res = aggregate(structs,defaultop,fieldops)
% we skip empty records
structs = structs(~cellfun('isempty',structs));

% for any struct array in the input, we first merge it recursively as if it were a cell array
for k=find(cellfun('length',structs)>1)
    structs{k} = hlp_aggregatestructs(structarray2cellarray(structs{k}),defaultop,fieldops); end

% at this point, we should have a cell array of structs
if ~isempty(structs)
    % we begin with the first struct
    res = structs{1};
    % and aggregate the remaining ones onto it
    for i=2:length(structs)
        si = structs{i};
        % proceeding field by field...
        for fn=fieldnames(si)'
            f = fn{1};
            % figure out which operation applies
            if isfield(fieldops,f)
                % a field-specific op applies
                if isstruct(fieldops.(f))
                    % ... which is itself a struct
                    fop = fieldops.(f);
                else
                    op = fieldops.(f);
                end
            else
                % the default op applies
                op = defaultop;
                fop = fieldops;
            end
            
            % now process the field
            if ~isfield(res,f)
                % field is not yet in the aggregate: just assign
                res.(f) = si.(f);
            else
                % need to aggregate it
                if isstruct(res.(f)) && isstruct(si.(f))
                    % both are a struct: recursively aggregate
                    res.(f) = aggregate({res.(f),si.(f)},op,fop);
                else
                    % they are not both structus
                    try
                        % try to apply the combiner op
                        res.(f) = op{1}(res.(f),si.(f));
                    catch
                        % didn't work: try to concatenate as fallback
                        try
                            res.(f) = [res.(f),si.(f)];
                        catch
                            % didn't work: try to assign as fallback (dropping previous field)
                            res.(f) = si.(f);
                        end
                    end
                end
            end
        end
    end
else
    % nothing to aggregate
    res = [];
end



function x = finalize(x,defaultop,fieldops)
% proceed field by field...
for fn=fieldnames(x)'
    f = fn{1};
    % figure out which operation applies
    if ~isempty(fieldops) && isfield(fieldops,f)
        % a field-specific op applies
        if isstruct(fieldops.(f))
            % ... which is itself a struct
            fop = fieldops.(f);
        else
            op = fieldops.(f);
        end
    else
        % the default op applies
        op = defaultop;
        fop = fieldops;
    end
    try
        % now apply the finalizer
        if isstruct(x.(f))
            % we have a sub-struct: recurse
            x.(f) = finalize(x.(f),op,fop);
        else
            % we have a regular element: apply finalizer
            x.(f) = op{2}(x.(f));
        end
    catch
        % for empty structs, x.(f) produces no output
    end
end



% translate string ops into actual ops, add the default finalizer if missing
function op = translateop(op)
if isstruct(op)
    % recurse
    op = structfun(@translateop,op,'UniformOutput',false);
else
    % remap strings
    if ischar(op)
        switch op
            case 'cat'
                op = @(a,b)[a b];
            case 'replace'
                op = @(a,b)b;
            case 'sum'
                op = @(a,b)a+b;
            case 'mean'
                op = {@(a,b)[a b], @(x)mean(x)};
            case 'median'
                op = {@(a,b)[a b], @(x)median(x)};
            case 'std'
                op = {@(a,b)[a b], @(x)std(x)};
            case 'random'
                op = {@(a,b)[a b], @(x) x(min(length(x),ceil(eps+rand(1)*length(x))))};
            case 'fillblanks'
                op = @(a,b)fastif(isempty(a),b,a);
            otherwise
                error('unsupported combiner op specified');
        end
    end
    % add finalizer if missing
    if ~iscell(op)
        op = {op,@(x)x}; end
end



% inefficiently turn a struct array into a cell array of structs
function res = structarray2cellarray(arg)
res = {};
for k=1:numel(arg)
    res = [res {arg(k)}]; end


% for the 'fillblanks' translate op
function val = fastif(cond,trueval,falseval)
if cond
    val = trueval;
else
    val = falseval;
end
