function hlpText = hlp_getFcnPreambleText(functionName)
% Return the 'preamble' text for a given function name.
% This is the initial part of the help text, and is terminated by a 
% horizontal dashed line at least 5 characters long. The preamble text 
% cannot contain any uncommented line breaks as these will be interpreted 
% as indicating the end of the block of help text.
%
% Trailing whitespace is removed from the preamble.
%
% Author: Tim Mullen June 2012, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


% make sure function name is valid
try
    utl_whichfile(functionName);
catch
    error('SIFT:hlp_getFcnPreambleText:badFcnName','Unable to find function %s.m',functionName);
end

% Next, read in the help text
if isdeployed
    hlpText = hlp_readHelpInDeployedMode(functionName);
else
    hlpText = help(functionName);
end

if isempty(hlpText)
    error('Error parsing %s.m: I could not find a block of help text for the function',functionName);
end

% parse the help text. We only want the text preceeding the first
% sufficiently long horizontal dashed line
endpoint = strfind(hlpText,'-----');

if isempty(endpoint)
    disp(hlpText)
    % error('Error parsing %s.m: I could not find the horizontal dashed line separating the preamble from remaining help text',functionName);
    warning('Error parsing %s.m: I could not find the horizontal dashed line separating the preamble from remaining help text',functionName);
    hlpText = '';
end

% return the preamble
hlpText = deblank(hlpText(1:endpoint-1));


