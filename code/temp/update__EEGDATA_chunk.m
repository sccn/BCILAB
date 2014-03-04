EEGDATA.buffer(:,1+mod(EEGDATA.smax:EEGDATA.smax+size(EEGDATA_chunk,2)-1,EEGDATA.buffer_len)) = EEGDATA_chunk;
EEGDATA.smax = EEGDATA.smax + size(EEGDATA_chunk,2);