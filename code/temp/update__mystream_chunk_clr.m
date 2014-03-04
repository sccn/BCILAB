mystream_range = 1+mod(mystream.smax:mystream.smax+size(mystream_chunk_clr,2)-1,mystream.buffer_len);
mystream.marker_pos(:,mystream_range) = 0;
mystream.buffer(:,mystream_range) = mystream_chunk_clr;
mystream.smax = mystream.smax + size(mystream_chunk_clr,2);