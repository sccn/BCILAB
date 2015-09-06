function com = hlp_cfg2comstr(config,mandatoryArgs,functionName)
% convert a config structure or cell array of NVPs to a command string that can be
% executed using eval
% 
% Author: Tim Mullen, 2012, SCCN/INC/UCSD

% Usage: to generate mandatoryArgs, can use inputname() to retrieve the
% name in the caller workspace of a variable. If the input is a literal
% (e.g. numeric) then inputname() will be '' and can just retrieve value
% from varargin

configOnly = nargin<2;
    
if ~configOnly && ~exist(functionName,'file')
    warning('function %s() does not exist',functionName)
end

if isstruct(config)
    % convert to cell array
    config = hlp_struct2vararginrec(config);
end

% clean up the cell array
config = hlp_cleanUp(config);

% convert cell array to evaluable string
argstr = hlp_vararg2str(config);


if configOnly
    % return only config NVPs as string
    com = argstr;
else
    % build command string of form 'functionName(mandatoryArgs,config);'
    com = sprintf('%s(%s,%s);',functionName,mandatoryArgs,argstr);
end



function cellarr = hlp_cleanUp(cellarr)

    % clean this depth level

    % remove any occurences of 'arg_direct' and paired (subsequent) value
    argpos = find(cellfun(@(x) isequal(x,'arg_direct'),cellarr));
    if ~isempty(argpos)
        argpos = [argpos argpos+1];
    end     
    cellarr(argpos) = [];

    % find any occurance of arg_selection and move it and subsequent value to
    % front of cell array
    argpos = find(cellfun(@(x) isequal(x,'arg_selection'),cellarr));
    if ~isempty(argpos)
        if isequal(cellarr{argpos+1},0)
            % arg_selection is 0 so replace entire cell array with []
            cellarr = [];
        	return;
        end
        cellarr = [cellarr(argpos:argpos+1) cellarr]; % ...copy values
        argpos  = argpos+2;
        cellarr(argpos:argpos+1) = []; % ...remove ('arg_selection', value)
    end

    % clean the next depth level
    for i=1:length(cellarr),
        if iscell(cellarr{i}),
            cellarr{i} = hlp_cleanUp(cellarr{i});
        end
    end
