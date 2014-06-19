background_range = 1+mod(background.smax:background.smax+size(background_chunk_clr,2)-1,background.buffer_len);
background.marker_pos(:,background_range) = 0;
background.buffer(:,background_range) = background_chunk_clr;
background.smax = background.smax + size(background_chunk_clr,2);