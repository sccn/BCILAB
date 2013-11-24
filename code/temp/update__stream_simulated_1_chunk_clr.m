stream_simulated_1_range = 1+mod(stream_simulated_1.smax:stream_simulated_1.smax+size(stream_simulated_1_chunk_clr,2)-1,stream_simulated_1.buffer_len);
stream_simulated_1.marker_pos(:,stream_simulated_1_range) = 0;
stream_simulated_1.buffer(:,stream_simulated_1_range) = stream_simulated_1_chunk_clr;
stream_simulated_1.smax = stream_simulated_1.smax + size(stream_simulated_1_chunk_clr,2);