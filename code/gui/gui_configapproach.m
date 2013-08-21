function [approach,action] = gui_configapproach(approach,dosave)
% bring up a modal configuration dialog for the given approach
% [Approach,Action] = gui_configapproach(Approach)
%
% In:
%   Approach : an approach; struct with fields 'paradigm' and 'parameters' (and optionally 'description' and 'name')
%              or cell array {paradigm, parameter1, parameter2, ...}
%
%   DoSave: Whether to bring up a save approach gui, after clicking okay (default: false)
%
% Out:
%   Result : a (re-)configured version of the Approach, or the unmodified input Approach (though possibly reformatted) if the user pressed 'Cancel'
%   Action : either 'OK' or 'Cancel'
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-10-25

if ~exist('approach','var') || isempty(approach)
    approach = gui_chooseapproach(); 
    if isempty(approach)
        return; end
end

if ~exist('dosave','var') || isempty(dosave)
    dosave = false; end

if iscell(approach)
    approach = struct('paradigm',approach{1}, 'parameters', {approach(2:end)}); end

calibrate_func = approach.paradigm;
if ischar(calibrate_func)
    instance = eval(calibrate_func); %#ok<NASGU>
    calibrate_func = eval('@instance.calibrate');
end
    
% bring up the GUI

result = arg_guidialog(calibrate_func,'params',approach.parameters,'title','BCILAB: Configure approach','Invoke',false);
if ~isempty(result)
    approach.parameters = {result};
    % get rid of hidden references...
    approach = utl_prune_handles(approach);
    if dosave
        gui_saveapproach(approach); end
    action = 'OK';
else
    action = 'Cancel';
end
