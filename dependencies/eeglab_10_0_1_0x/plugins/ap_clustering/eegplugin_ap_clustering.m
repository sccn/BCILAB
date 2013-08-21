% eegplugin_ap_clustering() - Add both Affinity Product pre-clustering AND Clustering commands to
% STUDY menu.
function vers = eegplugin_ap_clustering( fig, try_strings, catch_strings)
 
% title shown when loading into eeglab at the beginning
vers = 'Affinity_Product 1.0';

% create menu
plotmenu = findobj(fig, 'tag', 'study');

% make a sub-menu under study
submenu = uimenu( plotmenu, 'Label', 'Affinity Product clustering', 'position', 6, 'enable', 'on');

uimenu( submenu, 'label', 'Build preclustering matrices', 'callback', [ try_strings.no_check '[STUDYTMP ALLEEGTMP LASTCOM]= pop_apreclust(STUDY, ALLEEG);' catch_strings.update_study]);
uimenu( submenu, 'label', 'Cluster components', 'callback', [ try_strings.no_check '[STUDYTMP ALLEEGTMP LASTCOM] = pop_apcluster(STUDY, ALLEEG);' catch_strings.update_study]);
