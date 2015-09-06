function [signal] = flt_spect_features(varargin)
% Transform bandpower features
%

if ~exp_beginfun('filter') return; end

% 1:  [8 12]  ch1 
% 2:  [13 19] ch1 
% 3:  [2 4]   ch1 
% 4:  [4 7]   ch1
% 5:  [35 50] ch1
% 6:  [8 12]  ch2 
% 7:  [13 19] ch2 
% 8:  [2 4]   ch2 
% 9:  [4 7]   ch2
% 10: [35 50] ch2
% 11: [8 12]  ch3 
% 12: [13 19] ch3 
% 13: [2 4]   ch3 
% 14: [4 7]   ch3
% 15: [35 50] ch3
% 16: [8 12]  ch4 
% 17: [13 19] ch4 
% 18: [2 4]   ch4 
% 19: [4 7]   ch4
% 20: [35 50] ch4

rightChan = [1:10];
leftChan = [11:20];
alphaChan = [1 6 11 16];
betaChan = [2 7 12 17];
deltaChan = [3 8 13 18];
thetaChan = [4 9 14 19];
gammaChan = [5 10 15 20];

%
% eeg_chunk.data is a vector of bandpower features
% if single-channel or mean-channel:
%   dat = [band1 band2 ... bandN]
% else if k-channels
%   dat = [band1_ch1 ... bandN_ch1, band1_ch2 ... bandN_ch2, ... band1_chK ... bandN_chK]
%   

persistent state;

declare_properties('name', 'SpectralFeatures','follows','flt_fourier_bandpower','cannot_follow','set_makepos', 'independent_trials',true, 'independent_channels',false);

arg_define(varargin, ...
        arg_norep({'signal','Signal'}), ...
        arg_subtoggle('ztransform','on',@flt_zscore,'z-score spectral features'), ...
        arg({'featureOption1','Feature1'},'InterHemis_All_Valence',...
        {'none','InterHemis_Alpha_Valence','InterHemis_All_Valence',...
        'Alpha_Arousal','Theta_Valence','Lin2010_Valence',...
        'Lin2010_Arousal','Lin2013_Valence','Lin2013_Arousal'},'help1'), ...
        arg({'featureOption2','Feature2'},'Alpha_Arousal',...
        {'none','InterHemis_Alpha_Valence','InterHemis_All_Valence',...
        'Alpha_Arousal','Theta_Valence','Lin2010_Valence',...
        'Lin2010_Arousal','Lin2013_Valence','Lin2013_Arousal'},'help2'));
% 
% 
% {sprintf('Spectral Features\n'),[sprintf(['opt1: Interhemispheric Ratio over Alpha Band (Schmidt and Trainor)\n' ...
%                         'Left: Positive Valence, Right: Negative Valence\n\n']), ...
%          sprintf(['opt2: Interhemispheric Ratio over all bands (Altenmuller et al):\n' ...
%                          'Left: Positive Valence, Right: Negative Valence\n\n'])]}
%                      
%                      {sprintf('Spectral Features\n'),[sprintf(['opt1: Interhemispheric Ratio over Alpha Band (Schmidt and Trainor)\n' ...
%                         'Left: Positive Valence, Right: Negative Valence\n\n']), ...
%          sprintf(['opt2: Interhemispheric Ratio over all bands (Altenmuller et al):\n' ...
%                          'Left: Positive Valence, Right: Negative Valence\n\n'])]}
%                      
dat = signal.data;

switch featureOption1
    case 'none'
        datout = dat;
    case 'InterHemis_Alpha_Valence'
        %%% Interhemispheric Ratio over Alpha Band (Schmidt and Trainor): 
        %%% Left: Positive Valence, Right: Negative Valence
        datout(1) = sum(dat([1 5]))/(sum(dat([9 13])));
    case 'InterHemis_All_Valence'
        %%% Interhemispheric Ratio over all bands (Altenmuller et al): 
        %%% Left: Positive Valence, Right: Negative Valence
        datout(1) = (sum(dat(leftChan))/(sum(dat(rightChan))));
    case 'Alpha_Arousal'
        %%% Overall Alpha associated with arousal and fast music (Lin 2013)
        datout(1) = sum(dat(alphaChan))/4;
    case 'Theta_Valence'
        %%% Increase in Theta associated with positive valence 
        %%% (Sammler et al)
        datout(1) = sum(dat(thetaChan))/4;
    case 'Lin2010_Valence'
        %%% Increase in Theta, Decrease in Delta associated with positive 
        %%% valence (Lin 2010)
        datout(1) = (sum(dat(thetaChan))-sum(dat(deltaChan)));
    case 'Lin2010_Arousal'
        %%% Increase in Theta, Increase in Delta associated with increased 
        %%% arousal (Lin 2010)
        datout(1) = (sum(dat(thetaChan))+sum(dat(deltaChan)))/8;
    case 'Lin2013_Valence'
        %%% Increase in Gamma associated with major mode (Lin 2013)
        datout(1) = sum(dat(gammaChan))/4;
    case 'Lin2013_Arousal'
        %%% Increase in Alpha associated with fast music (Lin 2013) and
        %%% Increase in Theta associated with slow music
        datout(1) = (sum(dat(alphaChan))-sum(dat(thetaChan)));
end

switch featureOption2
    case 'none'
        
    case 'InterHemis_Alpha_Valence'
        %%% Interhemispheric Ratio over Alpha Band (Schmidt and Trainor): 
        %%% Left: Positive Valence, Right: Negative Valence
        datout(2) = sum(dat([1 5]))/(sum(dat([9 13])));
    case 'InterHemis_All_Valence'
        %%% Interhemispheric Ratio over all bands (Altenmuller et al): 
        %%% Left: Positive Valence, Right: Negative Valence
        datout(2) = (sum(dat(leftChan))/(sum(dat(rightChan))));
    case 'Alpha_Arousal'
        %%% Overall Alpha associated with arousal and fast music (Lin 2013)
        datout(2) = sum(dat(alphaChan))/4;
    case 'Theta_Valence'
        %%% Increase in Theta associated with positive valence 
        %%% (Sammler et al)
        datout(2) = sum(dat(thetaChan))/4;
    case 'Lin2010_Valence'
        %%% Increase in Theta, Decrease in Delta associated with positive 
        %%% valence (Lin 2010)
        datout(2) = (sum(dat(thetaChan))-sum(dat(deltaChan)));
    case 'Lin2010_Arousal'
        %%% Increase in Theta, Increase in Delta associated with increased 
        %%% arousal (Lin 2010)
        datout(2) = (sum(dat(thetaChan))+sum(dat(deltaChan)))/8;
    case 'Lin2013_Valence'
        %%% Increase in Gamma associated with major mode (Lin 2013)
        datout(2) = sum(dat(gammaChan))/4;
    case 'Lin2013_Arousal'
        %%% Increase in Alpha associated with fast music (Lin 2013) and
        %%% Increase in Theta associated with slow music
        datout(2) = (sum(dat(alphaChan))-sum(dat(thetaChan)));
end
        
if ztransform.arg_selection
    [signal state] = hlp_scope({'disable_expressions',true},@flt_zscore,'signal',signal,ztransform,'state',state);
end

signal.data = datout;

exp_endfun;