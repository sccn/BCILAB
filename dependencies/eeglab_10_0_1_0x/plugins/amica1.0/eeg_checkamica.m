% provides consistency between LLt and v, and between LLt_smooth and v
% stored under EEG.etc.amica
function EEG = eeg_checkamica(EEG)

if isfield(EEG.etc,'amica') && ~isempty(EEG.etc.amica)
    if isfield(EEG.etc.amica,'LLt') && ~isempty(EEG.etc.amica.LLt) 
        if ~isfield(EEG.etc.amica,'v') | ~isequal(size(EEG.etc.amica.LLt),size(EEG.etc.amica.v))
            EEG.etc.amica.v = LLt2v(EEG.etc.amica.LLt);
        end
    end
    if isfield(EEG.etc.amica,'LLt_smooth') && ~isempty(EEG.etc.amica.LLt_smooth)
        if ~isfield(EEG.etc.amica,'v_smooth') | ~isequal(size(EEG.etc.amica.LLt_smooth),size(EEG.etc.amica.v_smooth))
            EEG.etc.amica.v_smooth = LLt2v(EEG.etc.amica.LLt_smooth);
        end
    end
end