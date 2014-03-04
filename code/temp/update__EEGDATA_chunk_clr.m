EEGDATA_range = 1+mod(EEGDATA.smax:EEGDATA.smax+size(EEGDATA_chunk_clr,2)-1,EEGDATA.buffer_len);
EEGDATA.marker_pos(:,EEGDATA_range) = 0;
EEGDATA.buffer(:,EEGDATA_range) = EEGDATA_chunk_clr;
EEGDATA.smax = EEGDATA.smax + size(EEGDATA_chunk_clr,2);