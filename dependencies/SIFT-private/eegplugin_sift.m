% eegplugin_sift() - EEGLAB plugin for Source Information Flow Toolbox
%
% Usage:
%   >> eegplugin_sift(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   This plugins consist of the following Matlab files:
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Tim Mullen and Arnaud Delorme, SCCN, INC, UCSD

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2010 Tim Mullen
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
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

function vers = eegplugin_sift(fig, trystrs, catchstrs)

    siftroot = fileparts(which('StartSIFT.m'));
    
    fid = fopen([siftroot filesep 'resources' filesep 'version.txt']);
    version = fscanf(fid,'%s');
    fclose(fid);

    vers = num2str(version);
    
    if nargin < 3
        error('eegplugin_sift requires 3 arguments');
    end;
    
    % run SIFT startup routines
    ok=StartSIFT(false);
    
    if ~ok
        fprintf('SIFT initialization failed!\n');
        pause(2);
        return;
    end
    
    refreshMenuCB = 'hlp_refreshEEGLABMenus(EEG);';
    
    % find import data menu
    % ---------------------
    highlevelmenu = findobj(fig, 'tag', 'tools');
    
    % menu callbacks
    % --------------
    
    % simulation
    cmd = 'EEG = pop_sim_varmodel(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    SimVar_callback        = [finalcmd catchstrs.new_and_hist refreshMenuCB];
    
    cmd = 'EEG = pop_sim_eegdata(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    SimEEG_callback        = [finalcmd catchstrs.new_and_hist refreshMenuCB];
    
    % pre-processing
    cmd = 'EEG = pop_pre_prepData(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    PreProc_callback        = [finalcmd catchstrs.new_and_hist refreshMenuCB];
    
    % model fitting and validation
    cmd = 'EEG = pop_est_fitMVAR(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    FitModel_callback       = [finalcmd catchstrs.store_and_hist refreshMenuCB];
    
    cmd  = 'EEG = pop_est_selModelOrder(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    SelectModelOrder_callback   = [finalcmd catchstrs.store_and_hist refreshMenuCB];
    
    cmd  = 'EEG = pop_est_validateMVAR(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    ValidateModel_callback   = [finalcmd catchstrs.store_and_hist refreshMenuCB];
    
    % connectivity
    cmd = 'EEG = pop_est_mvarConnectivity(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    Connectivity_callback   = [finalcmd catchstrs.store_and_hist refreshMenuCB];

    % statistics
    cmd = 'EEG = pop_stat_surrogateGen(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    SurrogateDistrib_callback = [finalcmd catchstrs.store_and_hist refreshMenuCB];
    
    cmd = 'EEG = pop_stat_surrogateStats(EEG,0); CURRENTSET = CURRENTSET(1);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    SurrogateStats_callback = [finalcmd catchstrs.new_and_hist refreshMenuCB];
    
    cmd = 'EEG = pop_stat_analyticStats(EEG,0);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    AnalyticStat_callback   = [finalcmd catchstrs.store_and_hist refreshMenuCB];
    
    SimpleStat_callback     = 'warndlg2(''Coming Soon!'')';
    
    
    % visualization
    cmd = 'EEG = pop_vis_TimeFreqGrid(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    TFGrid_callback         = [finalcmd catchstrs.store_and_hist refreshMenuCB];
    
    CausalProjection_callback = 'warndlg2(''Coming Soon!'')';
    
    cmd = 'pop_vis_causalBrainMovie3D(EEG);';
    finalcmd = [ trystrs.no_check cmd ];
    finalcmd = [finalcmd 'LASTCOM = ''' cmd ''';' ];
    BranMovie_callback      = [finalcmd catchstrs.store_and_hist refreshMenuCB];
    
    % group analysis
    % ...
    
    % help
    AboutSIFT_callback  = 'gui_splashscreen;';
    SIFTManual_callback = @(hObject,eventdata) open([siftroot filesep 'documentation' filesep 'SIFT_manual.pdf']);
    SIFTWiki_callback   = 'web(''http://sccn.ucsd.edu/wiki/SIFT'',''-browser'');';
    
    
    % create menus
    % ------------
    userdata = 'startup:on;study:on';
    menu = uimenu( highlevelmenu, 'label', 'SIFT', 'separator', 'on' ,'userdata', userdata,'tag','siftroot','callback',refreshMenuCB);
    set(highlevelmenu,'enable','on','userdata','startup:on;study:on');
    
    % Simulation
    simmenu     = uimenu( menu, 'label', 'Simulation','userdata', userdata, 'enable','on');
        simmenu_LDS = uimenu( simmenu, 'label', 'Linear Dynamical System','userdata', userdata);
            uimenu( simmenu_LDS, 'label', 'Vector Autoregressive Process', 'callback', SimVar_callback,'userdata', userdata,'tag','sim');
        simmenu_EEG = uimenu( simmenu, 'label', 'EEG Simulator','callback', SimEEG_callback,'userdata', userdata,'tag','sim');
            
    % Preprocessing
    uimenu( menu, 'label', 'Pre-processing', 'callback', PreProc_callback,'userdata', userdata,'tag','cat');
    
    % Model fitting
    modelmenu   = uimenu( menu, 'label', 'Model fitting and validation','userdata', userdata);
        uimenu( modelmenu, 'label', 'Model Order Selection', 'callback', SelectModelOrder_callback ,'userdata', userdata,'tag','ic');
        uimenu( modelmenu, 'label', 'Fit AMVAR Model', 'callback',FitModel_callback,'userdata', userdata,'tag','model');
        uimenu( modelmenu, 'label', 'Validate model', 'callback', ValidateModel_callback ,'userdata', userdata,'tag','validation');
    
    % Connectivity
    connectmenu = uimenu( menu, 'label', 'Connectivity'  ,'callback',Connectivity_callback,'userdata', userdata,'tag','conn');
    
    % Statistics
    statmenu    = uimenu( menu, 'label', 'Statistics','userdata', userdata);
        uimenu( statmenu, 'label', 'Surrogate Distributions', 'callback', SurrogateDistrib_callback ,'enable','on','userdata', userdata,'tag','surogdist');
        uimenu( statmenu, 'label', 'Surrogate Statistics', 'callback', SurrogateStats_callback ,'enable','on','userdata', userdata,'tag','stats');
        uimenu( statmenu, 'label', 'Analytic Statistics (Beta)', 'callback', AnalyticStat_callback,'separator', 'on', 'enable','on','userdata', userdata);
        %uimenu( statmenu, 'label', 'Simple Statistics', 'callback', SimpleStat_callback, 'separator', 'on' ,'enable','off');
    
    % Visualization
    vismenu     = uimenu( menu, 'label', 'Visualization' ,'userdata', userdata);
        uimenu( vismenu , 'label', 'Time-Frequency Grid', 'callback', TFGrid_callback ,'userdata', userdata);
        uimenu( vismenu , 'label', 'BrainMovie3D', 'callback', BranMovie_callback ,'userdata', userdata);
        %uimenu( vismenu , 'label', 'Causal Projection', 'callback', CausalProjection_callback, 'enable','off' );
    
    % Group Analysis
    grpmenu     = uimenu( menu, 'label', 'Group Analysis' ,'userdata', userdata);
        grpmenu_mpa = uimenu( grpmenu , 'label', 'Causal Projection (MPT)','userdata', 'study:on');
        set(grpmenu,'callback',{@create_mpt_menus,grpmenu_mpa});
        
    % Help
    helpmenu    = uimenu( menu, 'label', 'Help','userdata', userdata, 'enable','on');
        aboutmenu   = uimenu( helpmenu, 'label', 'About SIFT' ,'callback',AboutSIFT_callback,'userdata', userdata);
        manualmenu  = uimenu( helpmenu, 'label', 'SIFT Manual' ,'callback',SIFTManual_callback,'userdata', userdata , 'enable','on');
        wikimenu    = uimenu( helpmenu, 'label', 'SIFT Wiki' ,'callback',SIFTWiki_callback,'userdata', userdata ,'enable','on');

 
    % refresh the checkmark state of EEGLAB menus
    eeg_global;
    hlp_refreshEEGLABMenus(EEG);
    
end

% Helper Functions (callbacks)
% ------------------------------------------------------------------------------------------------------
function create_mpt_menus(h,ev,grpmenu_mpa)
    if ~isempty(get(grpmenu_mpa,'children')) return; end

     % global variables that are used by measure projection
    global visualize_menu_label
    global visualize_by_measure_menu_label
    global visualize_by_domain_menu_label
    global visualize_by_significance_menu_label

    if exist('eegplugin_mproject','file')
        set(grpmenu_mpa,'enable','on');

        visualize_menu_label = 'Show volume';
        visualize_by_measure_menu_label = 'Show colored by Measure';
        visualize_by_domain_menu_label = 'Show colored by Domain';
        visualize_by_significance_menu_label =  'Show colored by Significance';

        % create measure projection submenus
        create_mpt_submenu(grpmenu_mpa, 'sift');  
        uimenu(grpmenu_mpa, 'Label', 'Options', 'Callback', @show_gui_options, 'separator','on','userdata', 'study:on');
        uimenu(grpmenu_mpa, 'Label', 'About', 'Callback', @show_about, 'separator','on', 'userdata', 'study:on');

        update_domain_menu([],[],grpmenu_mpa);
        set(grpmenu_mpa,'callback',@(h,tmp) update_domain_menu(h,tmp,grpmenu_mpa));
    else
        set(grpmenu_mpa,'enable','off');
    end
end
        
function update_domain_menu(h,ev,mprojection_SIFTmenu)
    
    maxDomainNo = 30;

    try
        studyHasSiftDomain = evalin('base', '~isempty(STUDY.measureProjection.sift.projection.domain)');
    catch
        studyHasSiftDomain = false;
    end;

    mainfig = findobj('tag','EEGLAB');
    if ~isempty(mainfig)
        siftDomainMenu = findobj('parent', mprojection_SIFTmenu , 'type', 'uimenu', 'label','Domains');
        set(siftDomainMenu,'enable',fastif(studyHasSiftDomain,'on','off'));

        if studyHasSiftDomain
            nbSiftDomains = evalin('base','size(STUDY.measureProjection.sift.projection.domain,2)');
            for i = 1:maxDomainNo
                if i <= nbSiftDomains
                    set(findobj('parent',siftDomainMenu,'Label',['Domain ' num2str(i)]),'visible','on');
                else
                    set(findobj('parent',siftDomainMenu,'Label',['Domain ' num2str(i)]),'visible','off');

                end
            end
        end

    else
        evalin('base','eeglab redraw');
        update_domain_menu;
    end
end