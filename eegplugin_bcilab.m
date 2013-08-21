function vers = eegplugin_bcilab(fig,trystrs,catchstrs)
% BCILAB plugin function for EEGLAB
%
% This file will only be used if you place the BCILAB folder in plugins folder of an EEGLAB installation.
%
% Notes:
%  The file bcilab_config.m may contain definitions of env_startup parameters; these will be automatically applied when the plugin is loaded.
%
%  BCILAB can also be loaded as a standalone toolbox with EEGLAB as a dependency, by removing the BCILAB folder from the EEGLAB/plugins folder, and putting
%  the entire EEGLAB folder as-is into the dependencies folder of BCILAB.
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-10-30

% invoke the configuration
bcilab_config

% collect the variables as defined in bcilab_config in a struct
infos = whos();
for n = {infos.name}
    settings.(n{1}) = eval(n{1}); end

% load the toolbox
cd([fileparts(mfilename('fullpath')) filesep 'code' filesep 'environment']); 
env_startup(settings);

% return the version
vers = env_version;
