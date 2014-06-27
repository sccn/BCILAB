try
    EEGDATA.marker_buffer(:,1+mod(EEGDATA.mmax:EEGDATA.mmax+size(EEGDATA_marker_chunk,2)-1,EEGDATA.marker_buffer_len)) = EEGDATA_marker_chunk;
catch e
    if any(strcmp(e.identifier,{'MATLAB:heterogeneousStrucAssignment','MATLAB:heterogenousStrucAssignment'}))
        EEGDATA.marker_buffer = EEGDATA_marker_chunk(ones(1,EEGDATA.marker_buffer_len));
        EEGDATA.marker_buffer(:,1+mod(EEGDATA.mmax:EEGDATA.mmax+size(EEGDATA_marker_chunk,2)-1,EEGDATA.marker_buffer_len)) = EEGDATA_marker_chunk;
    else
        rethrow(e);
    end
end