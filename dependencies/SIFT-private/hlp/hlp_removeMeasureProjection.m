function hlp_removeMeasureProjection

% temporary hack to remove measure projection toolbox from path 
% (MPT interferes with SIFT)
p = which('eegplugin_sift.m');
p = p(1:findstr(p,'eegplugin_sift.m')-1);
curdir = pwd;

cd([p '../']);
mp_paths = genpath(fullfile(pwd,'measure_projection'));
warn = warning;
warning off all
rmpath(mp_paths);
warning(warn);
cd(curdir);
% --