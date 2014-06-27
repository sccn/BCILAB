laststream_range = 1+mod(laststream.smax:laststream.smax+size(laststream_chunk_clr,2)-1,laststream.buffer_len);
laststream.marker_pos(:,laststream_range) = 0;
laststream.buffer(:,laststream_range) = laststream_chunk_clr;
laststream.smax = laststream.smax + size(laststream_chunk_clr,2);