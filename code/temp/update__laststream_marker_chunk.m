try
    laststream.marker_buffer(:,1+mod(laststream.mmax:laststream.mmax+size(laststream_marker_chunk,2)-1,laststream.marker_buffer_len)) = laststream_marker_chunk;
catch e
    if any(strcmp(e.identifier,{'MATLAB:heterogeneousStrucAssignment','MATLAB:heterogenousStrucAssignment'}))
        laststream.marker_buffer = laststream_marker_chunk(ones(1,laststream.marker_buffer_len));
        laststream.marker_buffer(:,1+mod(laststream.mmax:laststream.mmax+size(laststream_marker_chunk,2)-1,laststream.marker_buffer_len)) = laststream_marker_chunk;
    else
        rethrow(e);
    end
end