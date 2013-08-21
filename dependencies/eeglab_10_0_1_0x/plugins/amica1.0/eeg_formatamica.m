function EEG = eeg_formatamica(EEG)
if isfield(EEG.etc.amica,'LLt') && isfield(EEG.etc.amica,'LLt_smooth')
    if size(EEG.data,2) == size(EEG.etc.amica.LLt,2) && size(EEG.data,3) == size(EEG.etc.amica.LLt,3)
        if size(EEG.data,2) == size(EEG.etc.amica.LLt_smooth,2) && size(EEG.data,3) == size(EEG.etc.amica.LLt_smooth,3)
            
            EEG.data = [EEG.data; EEG.etc.amica.LLt; EEG.etc.amica.LLt_smooth];
            EEG.etc.amica.prob_added = 1;
            num_models = size(EEG.etc.amica.LLt,1);
            EEG.nbchan = EEG.nbchan + 2*num_models;
            if ~isempty(EEG.chanlocs)
                channames = {};
                for k=1:2*num_models
                    channames{end+1} = ['prob' num2str(k)]; end
                [EEG.chanlocs(end+1 : end+2*num_models).labels] = channames{:};
            end
        end
    end
end
