function EEG = fixegilocs(EEG, fileloc)

alllocs = readlocs(fileloc);

EEG.chaninfo.ndchanlocs = struct([]);

if (EEG.nbchan == 128 || EEG.nbchan == 256) && length(alllocs) == EEG.nbchan+4
    ndchanlocs = [1 2 3 length(alllocs)]; %3 fiducial and 1 reference channel
    
    alllocs(1).type = 'FID';
    alllocs(2).type = 'FID';
    alllocs(3).type = 'FID';
    alllocs(end).type = 'REF';
    
    EEG.chanlocs = alllocs(setdiff(1:length(alllocs),ndchanlocs));
    EEG.chaninfo.ndchanlocs = alllocs(ndchanlocs);
    EEG.chaninfo.filename = fileloc;
elseif EEG.nbchan == 7 && length(alllocs) == 7
    EEG.chanlocs = alllocs;
else
    error('Mismatch in channel counts: %d in data, %d in %s.', EEG.nbchan, length(alllocs), fileloc);
end

EEG = eeg_checkset(EEG);

fnlist = fieldnames(EEG.chaninfo.ndchanlocs);
for fn = 1:length(fnlist)
    if ~isfield(EEG.chanlocs,fnlist{fn})
        EEG.chaninfo.ndchanlocs = rmfield(EEG.chaninfo.ndchanlocs,fnlist{fn});
    end
end