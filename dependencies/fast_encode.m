function str = fast_encode(bytes)
% encode a uint8 vector into a string
bytes = uint8(32) + [bitand(bytes,uint8(15)) bitshift(bitand(bytes,uint8(15*16)),-4)];
str = char(bytes(:)');
