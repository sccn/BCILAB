function expose_handles(varargin)
% Helper function to make nested function handles accessible to arg_report.
% expose_handles(KnownHandles..., Varargin...)
%
% In some cases outside access to a function Func's nested or scoped functions is necessary, for
% example when a such handles are returned by Func and later serialized and deserialized using
% hlp_deserialize. This access is provided by adding an expose_handles line at the top of Func's
% body and passing all of Func's arguments into expose_handles.
%
% In:
%   KnownHandles... : known function handles to expose, if any
%   Varargin... : the varargin of the calling function
%
% Reports:
%   A cell array of looked-up function handles, or []'s for handles that could not be looked up
%   successfully.
%
% Examples:
%   function myfunc(a,b,varargin)
%   % (Your documentation...)
%
%   % expose nested/scoped function handles to the outside
%   expose_handles(a,b,varargin{:});
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-10-19

if length(varargin)>1 && ischar(varargin{end-1}) && strcmp(varargin{end-1},'__arg_report__') && ischar(varargin{end}) && any(strcmp(varargin{end},{'handle','handles'}))
    varargin = varargin(1:end-2);
    queries = varargin(cellfun('isclass',varargin,'char'));    
    if ~isempty(queries)
        handles = varargin(cellfun('isclass',varargin,'function_handle'));
        handle_strs = cellfun(@char,handles,'UniformOutput',false);
    end
    % try to resolve each function name in the caller scope
    for f=length(queries):-1:1
        if ~isempty(handle_strs) && any(strcmp(queries{f},handle_strs))
            funcs{f} = handles{strcmp(queries{f},handle_strs)}; 
        else
            % note: this code path seems to have stopped working in recent MATLABs
            funcs{f} = evalin('caller',['@' varargin{f}]); 
        end
        % check if lookup was successful
        if isempty(getfield(functions(funcs{f}),'file'))
            funcs{f} = []; end
    end
    % yield the report
    arg_issuereport(funcs);
end
