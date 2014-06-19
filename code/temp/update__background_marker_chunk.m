try
    background.marker_buffer(:,1+mod(background.mmax:background.mmax+size(background_marker_chunk,2)-1,background.marker_buffer_len)) = background_marker_chunk;
catch e
    if any(strcmp(e.identifier,{'MATLAB:heterogeneousStrucAssignment','MATLAB:heterogenousStrucAssignment'}))
        background.marker_buffer = background_marker_chunk(ones(1,background.marker_buffer_len));
        background.marker_buffer(:,1+mod(background.mmax:background.mmax+size(background_marker_chunk,2)-1,background.marker_buffer_len)) = background_marker_chunk;
    else
        rethrow(e);
    end
end