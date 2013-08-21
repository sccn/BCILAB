function success = utl_update_config(varargin)
% Update a set of values in the configuration file and handle permissions errors, etc.
% utl_update_config(NVPs...)
%
% In:
%   Operation : operation to perform on the config file; should be 'set'
%
%   NVPs... : list of name-value pairs, where each name denotes a config variables and the subsequent
%             value is the string expression that should be written into the config file. It is 
%             generally a good idea to use hlp_tostring() to turn a data structure into such a string
%             representation.
%
% Out:
%   Success : whether the update was successful
%
% Notes:
%   This function will try to update the currently active config file, which is by default the one
%   in the toolbox directory. If that file is read-only, the user will be asked whether he/she wants
%   to create their own personal copy of the config file (which is in the user's home directory, in a 
%   sub-folder named .bcilab.
%
% See also:
%   hlp_config
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-12-16

global tracking;
success = false;
try
    % apply updates
    hlp_config(tracking.configscript,varargin{:});
catch e
    if strcmp(e.identifier,'hlp_config:permissions_error')
        % no permission to update config file: ask if a local copy should be created
        resp = questdlg2('The configuration file that you use is read only. Would you like to create a local copy in your home directory?','Permissions','Yes','No','Cancel','Yes');
        if strcmp(resp,'Yes')
            try 
                % make writable copy
                new_config = env_translatepath('home:/.bcilab/bcilab_config.m');
                if ~exist(fileparts(new_config),'dir')
                    mkdir(env_translatepath('home:/'),'.bcilab'); end
                copyfile(tracking.configscript,new_config);                
                fileattrib(new_config,'+w');
                % make it the current config file
                tracking.configscript = new_config;
                % re-apply changes
                hlp_config(tracking.configscript,varargin{:});
            catch
                warndlg2(['Cannot create the file "' new_config '".'],'Notification');
                return;
            end
        end
    else
        warndlg2('Your new values cannot be applied due to some error.','Notification');
        env_handleerror(e);
        return;
    end
end

success = true;