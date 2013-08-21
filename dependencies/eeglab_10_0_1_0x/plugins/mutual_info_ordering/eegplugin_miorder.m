function eegplugin_miorder(fig, try_strings, catch_strings)


    % ---------------------
    toolsmenu = findobj(fig, 'tag', 'tools');
    submenu = uimenu(toolsmenu, 'label', 'Pairwise Mutual Information of ICs (beta)');
    uimenu( submenu, 'Label', 'Show pairwise MI', 'CallBack', ...
        [try_strings.check_ica 'pop_mishow(EEG);' catch_strings.store_and_hist '[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);' ]);
    uimenu( submenu, 'Label', 'Group ICs by MI' , 'CallBack', ...
        [try_strings.check_ica '[MInew ord] = pop_miorder2(EEG);' catch_strings.store_and_hist '[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);' ]);
    
    

    