function expose_handles(varargin)
% Helper function to make nested function handles accessible to arg_report.
% expose_handles(Varargin...)
%
% In some cases outside access to a function Func's nested or scoped functions is necessary, for
% example when a such handles are returned by Func and later serialized and deserialized using
% hlp_deserialize. This access is provided by adding an expose_handles line at the top of Func's
% body and passing all of Func's arguments into expose_handles.
%
% In:
%   Varargin : the varargin of the calling function
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

if length(varargin)>1 && strcmp(varargin{end-1},'__arg_report__') && strcmp(varargin{end},'handles')
    varargin = varargin(1:end-2);
    if ~iscellstr(varargin)
        error('The arguments passed for a handle report must denote function names.'); end
    % try to resolve each function name in the caller scope
    for f=length(varargin):-1:1
        funcs{f} = evalin('caller',['@' varargin{f}]); 
        % check if lookup was successful
        if isempty(getfield(functions(funcs{f}),'file'))
            funcs{f} = []; end
    end
    % yield the report
    arg_issuereport(funcs);
end
