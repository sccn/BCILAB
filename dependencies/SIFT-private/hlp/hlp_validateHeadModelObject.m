
function hmObj = hlp_validateHeadModelObject(hmObj)
% helper function to check the head model object, and
% load from base workspace if necessary
% Author: Tim Mullen, SCCN/INC/UCSD 2013

if isempty(hmObj)
    hmObj = 'hmObj'; % default name of variable in base wkspace
end
if ischar(hmObj)
    % load variable from base workspace
    try
       hmObj = evalin('base',hmObj);
    catch err1
        % try to open it as a file names
        try
            if ~exist(hmObj,'file'), hmObj = env_translatepath(hmObj); end
            hmObj = hlp_microcache('headmodels',@headModel.loadFromFile,hmObj);
        catch err2
            disp_once('BCILAB:HeadmodelMissing',['Note: The head model named "' hmObj '" was neither found in the base workspace nor could it be loaded as a file.']);
            hmObj = '';
            return;
        end
    end
end

% check that the object is a valid head model object
if ~isobject(hmObj)
    warning('You did not supply a valid head model object');
    hmObj = '';
    return;
end
if ~isprop(hmObj,'leadFieldFile') || isempty(hmObj.leadFieldFile)
    warning('Head model object does not appear to be valid. The lead field matrix does not exist.');
    hmObj = '';
    return;
end