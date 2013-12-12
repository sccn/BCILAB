try
    stream_simulated_1.marker_buffer(:,1+mod(stream_simulated_1.mmax:stream_simulated_1.mmax+size(stream_simulated_1_marker_chunk,2)-1,stream_simulated_1.marker_buffer_len)) = stream_simulated_1_marker_chunk;
catch e
    if any(strcmp(e.identifier,{'MATLAB:heterogeneousStrucAssignment','MATLAB:heterogenousStrucAssignment'}))
        stream_simulated_1.marker_buffer = stream_simulated_1_marker_chunk(ones(1,stream_simulated_1.marker_buffer_len));
        stream_simulated_1.marker_buffer(:,1+mod(stream_simulated_1.mmax:stream_simulated_1.mmax+size(stream_simulated_1_marker_chunk,2)-1,stream_simulated_1.marker_buffer_len)) = stream_simulated_1_marker_chunk;
    else
        rethrow(e);
    end
end