function result = verLessThan(toolboxstr, verstr)
%verLessThan Compare version of toolbox to specified version string.
%   verLessThan(TOOLBOX_DIR, VERSION) returns true if the version of
%   the toolbox specified by the string TOOLBOX_DIR is older than the
%   version specified by the string VERSION, and false otherwise. 
%   VERSION must be a string in the form 'major[.minor[.revision]]', 
%   such as '7', '7.1', or '7.0.1'. If TOOLBOX_DIR cannot be found
%   on MATLAB's search path, an error is generated.
%
%   Examples:
%       if verLessThan('images', '4.1')
%           error('Image Processing Toolbox 4.1 or higher is required.');
%       end
%
%       if verLessThan('matlab', '7.0.1')
%           % Put code to run under MATLAB older than MATLAB 7.0.1 here
%       else
%           % Put code to run under MATLAB 7.0.1 and newer here
%       end
%
%   See also MATLABPATH, VER.

%   Copyright 2006-2007 The MathWorks, Inc.
%   $Revision: 1.1.4.1 $  $Date: 2007/02/23 13:28:34 $
    
if nargin < 2
    errstr = 'Not enough input arguments.';
    if errorSupportsIdentifiers
        error('MATLAB:nargchk:notEnoughInputs', errstr)
    else
        error(errstr)
    end
end

if ~ischar(toolboxstr) | ~ischar(verstr)
    errstr = 'Inputs must be strings.';
    if errorSupportsIdentifiers
        error('MATLAB:verLessThan:invalidInput', errstr)
    else
        error(errstr)
    end
end

toolboxver = ver(toolboxstr);
if isempty(toolboxver)
    errformat = 'Toolbox ''%s'' not found.';
    if errorSupportsIdentifiers
        error('MATLAB:verLessThan:missingToolbox', errformat, toolboxstr)
    else
        error(sprintf(errformat, toolboxstr))
    end
end

toolboxParts = getParts(toolboxver(1).Version);
verParts = getParts(verstr);

result = (sign(toolboxParts - verParts) * [1; .1; .01]) < 0;

function parts = getParts(V)
    parts = sscanf(V, '%d.%d.%d')';
    if length(parts) < 3
       parts(3) = 0; % zero-fills to 3 elements
    end

function tf = errorSupportsIdentifiers
    % Determine, using code that runs on MATLAB 6.0 or later, if
    % error identifiers should be used when calling error().
    tf = 1;
    eval('lasterr('''','''');','tf = 0;');
