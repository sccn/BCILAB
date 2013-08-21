function pop_mishow(EEG)

if nargin < 1
    help pop_mishow;
    return;
end;

% check whether ica data is availible
% ---------------------------------------------------
EEG.icaact = eeg_getica(EEG);
if ~isempty( EEG.icaact )
    [dummy MI] = minfo(EEG.icaact);
    MI = (MI + MI') / 2;
    MI = MI - diag(diag(MI));
    
    figure; imagesc(MI);
else
    error('You must run ICA first');
end;