function EEG = eeg_reformatamica(EEG)
if isfield(EEG.etc.amica,'LLt')
    num_models = size(EEG.etc.amica.LLt,1);
    EEG.etc.amica.LLt_smooth = EEG.data(end-num_models+1:end,:,:);
    EEG.data(end-num_models+1:end,:,:) = [];
    EEG.etc.amica.LLt = EEG.data(end-num_models+1:end,:,:);
    EEG.data(end-num_models+1:end,:,:) = [];
    
    EEG.nbchan = EEG.nbchan - 2*num_models;
    if ~isempty(EEG.chanlocs)
        for i = 1:2*num_models;
            EEG.chanlocs(end) = [];
        end
    end
    EEG.etc.amica = rmfield(EEG.etc.amica,'prob_added');
end