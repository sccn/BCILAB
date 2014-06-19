background.buffer(:,1+mod(background.smax:background.smax+size(background_chunk,2)-1,background.buffer_len)) = background_chunk;
background.smax = background.smax + size(background_chunk,2);