function bytes = fast_decode(str)
bytes = uint8(str)-uint8(32);
bytes = (bytes(1:end/2) + bitshift(bytes((end/2+1):end),uint8(4)))';
