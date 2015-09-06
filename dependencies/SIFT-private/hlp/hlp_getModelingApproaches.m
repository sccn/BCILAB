function arglist = hlp_getModelingApproaches(varargin)
% return a cell array (arglist) defining the available modeling approaches
% each entry of the cell array is a cell array of the form
% {'Approach Name' @function_handle}
% where @function_handle is the handle to the function implementing the
% approach
%
% hlp_getModelingApproaches('defNameOnly') returns only the name of the 
%  default modeling approach (as a string)
%
% hlp_getModelingApproaches('mfileNameOnly',approachName) where approachName 
%   is a string with the name of a valid approach, returns the m-file name 
%   of the function implementing the specified approach
%
% Author: Tim Mullen, 2012 SCCN/INC/UCSD

defNameOnly = false;
approachName = '';

if nargin==1 && strcmp(varargin{1},'defNameOnly')
    defNameOnly = true;
end

if nargin==2 && ismember_bc('mfileNameOnly',varargin);
    approachName = varargin{2};
    if ~ischar(approachName)
        error('hlp_getModelingApproaches:badInput','Bad argument pair for ''mfileNameOnly'''); 
    end
end
    
% construct the arglist
arglist = { ...
            {'Segmentation VAR' @est_fitMVAR}, ...   % {'Kalman Filtering' {arg('arg_dummy',[],[],'replace this with @est_fitMVARKalman')}}, ...
            };
      

% determine the output
if defNameOnly
    arglist = arglist{1}{1};
elseif ~isempty(approachName)
    % get the function handle matching the desired approach name
    tmp = cellfun(@(x) fastif(isequal(x{1},approachName),x{2},''),arglist,'UniformOutput',false);
    tmp(cellfun(@isempty,tmp))=[];
    if isempty(tmp)
        error('hlp_getModelingApproachMfileName:badApproachName','Unknown modeling approach ''%s''',approachName);
    else
        % return function name as a string
        arglist = char(tmp{1});
    end
end
    