function result = bs_read(h)
% Read a block of new data from a BioSemi connection
% [Block] = bs_read(Handle)
%
% In:
%   Handle : handle to a BioSemi connection (previously opened via bs_open)
%
% Out:
%   Block : a block of new data, [#Channels x #Samples]
%           for EEG, the scale is in uV
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-12-09

% we keep track of the last buffer pointer in this variable
persistent last_ptr;
if length(last_ptr) < h.uid
    last_ptr(h.uid) = h.last_ptr; end

[raw,last_ptr(h.uid)] = bsb_read(h.hDLL,h.hConn,h.hREAD_POINTER,h.pBuffer,int32(h.nbchan),int32(last_ptr(h.uid)));

result = (single(raw))/8192;
