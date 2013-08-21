function str = eeglaboptionsfolder
% Return the folder in which EEGLAB should look for options

try
    % used with BCILAB
    if exist('env_translatepath','file')
        % during startup, if EEGLAB is loaded as dependency
        p = env_translatepath([pwd '/../../resources/eeg_optionsbackup.txt']);
        if exist(p,'file')
            str = fileparts(p);
            return;
        end
        % during operation, if EEGLAB is a dependency
        p = env_translatepath('bcilab:/resources/eeg_optionsbackup.txt');
        if exist(p,'file')
            str = fileparts(p);
            return;
        end
        % if BCILAB is a plugin to EEGLAB
        str = utl_whichfile('eeglab');
        if exist(p,'file')
            str = fileparts(p);
            return;
        end
    else
        % EEGLAB is the master app
        str = ctfroot;
        inds = find(str == filesep);
        str = str(1:inds(end)-1);
    end
catch
    % fall back original codepath
    str = ctfroot;
    inds = find(str == filesep);
    str = str(1:inds(end)-1);    
end
